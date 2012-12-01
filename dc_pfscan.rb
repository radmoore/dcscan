#!/usr/bin/ruby1.9.1

require 'thread'
require 'yaml'
require 'fileutils'

# dc_pfscan.rb
#
# A simple threaded wrapper around the pfam_scan utility for scanning
# against very large sequence files. Uses a divide and conquer approach.
# First splits large fasta into chunks of a defined max. size, then performs 
# the scan on these chunks in seperate threads using a job queue.
# Result files of each chunk are merged when the queue is empty.
# Configuration is done with a YAML file.
#
# More info on pfam_scan:
# ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/
# 
# This is part of the RADS annotation pipeline
#
# NOTE:
# The consequence of running seperate threads is that the order of
# the output may be jumbled (that is, as the runtime of pfam_scan
# will depend on sequence length & composition, number of hits etc
# such that chunks that were started later may be finished earlier.)
# That is to say that the *order* of hits obtained when running dc_pfscan
# and pfam_scan may be different.


# Mixin to count number of digits
# (for clean output)
class Integer
  def digit_no
    self.to_s.size
  end
end

# READ YAML
def read_config
  raw_config = File.read(ARGV[0])
  config = YAML.load(raw_config)
  if ( config[:files][:fasta].nil? )
    STDERR.puts "At least a fasta sequence must be provided (see #{ARGV[0]})"
    exit(-1)
  end
  return config
end

# SPLIT FASTA into CHUNKS
def split_fasta(config)
  fadir = "#{config[:files][:wdir]}/#{config[:files][:fadir]}"
  chunk_size = config[:jobcontrol][:chunk_size]
  width = chunk_size.digit_no
  if config[:files][:wipe] 
    if File.exists?(fadir) 
       FileUtils.rm_rf fadir
    end
  end
  Dir.mkdir(fadir) unless File.exists?(fadir)
  fasta = config[:files][:fasta]
  outfiles = []
  fasta_basename = File.basename(fasta, ".*") 
  in_file_no = 0
  out_file = File.new("#{fadir}/#{fasta_basename}-#{outfiles.size}.fa", "w")
  puts "\nA. PERFORMING SPLIT"
  puts "="*25
  print "\rSplitting [part %8d] %#{width}d... " % [outfiles.size, in_file_no]
  IO.foreach(fasta) do |line|
    if (line[0] == '>')
      if (in_file_no >= chunk_size)
        out_file.puts ""
        out_file.close()
        outfiles << out_file
        out_file = File.new("#{fadir}/#{fasta_basename}-#{outfiles.size}.fa", "w")
        puts "\rSplitting [part %8d] %#{width}d... done." % [outfiles.size, in_file_no]
        in_file_no = 0
      end
      in_file_no += 1
      print "\rSplitting [part %8d] %#{width}d... " % [outfiles.size, in_file_no]
    end
    out_file.puts line
  end
  out_file.close()
  outfiles << out_file
  print "\rSplitting [part %8d] %#{width}d... " % [outfiles.size, in_file_no]
  puts "done."
  return outfiles
end

# PERFORM SCAN ON CHUNKS
def run_pfamscan(infiles, config)
  puts "\nB. PERFORMING SCAN"
  puts "="*25
  pfsdir = "#{config[:files][:wdir]}/#{config[:files][:pfsdir]}"
  if config[:files][:wipe] 
    if File.exists?(pfsdir) 
       FileUtils.rm_rf pfsdir
    end 
  end 
  Dir.mkdir(pfsdir) unless File.exists?(pfsdir)
  threads = []
  outfiles = []
  exec_queue = Queue.new
  infiles.each do |f|
    base = File.basename(f.path, ".fa")
    outfile = "#{pfsdir}/#{base}.pfsout"
    exec_queue << "#{config[:pfamscan][:bin]} -fasta #{f.path} -dir #{config[:files][:dir]} #{(config[:pfamscan][:clans]) ? '-clan_overlap' : ''} -cpu #{config[:pfamscan][:cpu]} -outfile #{outfile}"
    outfiles << outfile
  end

  total_runs = outfiles.size
  with = total_runs.digit_no
  current_run = 0

  config[:jobcontrol][:max_thread].times do 
    threads << Thread.new do 
      until exec_queue.empty?
        scan_cmd = exec_queue.pop(true) rescue nil
        width = total_runs.digit_no
        puts "Scanning chunk %#{width}d of %d\n" % [current_run+=1, total_runs]
        `#{scan_cmd}` if not scan_cmd.nil?
        # not sure about this one
        if ($?.exitstatus != 0)
          STDERR.puts "There was a problem running #{config[:pfamscan][:bin]}."
          exit(-1)
        end
      end 
    end
  end
  threads.each{|thread| thread.join}
  return outfiles
  
end

# MERGE SCAN RESULTS
def merge_results(infiles, config)
  puts "\nC. MERGING RESULTS"
  puts "="*25
  if (config[:files][:resultfile].nil?)
    final_results = File.basename(config[:files][:fasta], ".*")
    outfile = File.new("#{final_results}.pfsout", "w")
  else
    outfile = File.new(config[:files][:resultfile], "w")
  end
  first = true
  width = infiles.size.digit_no
  current_merge = 0
  STDERR.print "Merging result files... "
  infiles.each do |f|
    IO.foreach(f) do |line|
      next if /^$/.match(line)
      next if line[0] == '#' && (not first)
      outfile.puts line
    end
    print "\rMerging result %#{width}d of %d... " % [current_merge+=1, infiles.size]
    first = false
  end
  puts "\rMerging result %#{width}d of %d... done." % [current_merge, infiles.size]
  puts ""
  outfile.close()
  return outfile.path
end

# CREATE XDOM FILE FROM SCAN RESULTS
def create_xdom(result_file, config)
  puts "\nD. CREATING XDOM FILE"
  puts "="*25
  if (config[:parsepfs][:filename].nil?)
    xdom_file = File.basename(config[:files][:fasta], ".*")
    outfile = File.new("#{xdom_file}.xdom", "w")
  else
    outfile = File.new(config[:parsepfs][:filename], "w")
  end

  dom_start_field = (config[:parsepfs][:envelope]) ? 2 : 3
  dom_end_field = (config[:parsepfs][:envelope]) ? 3 : 4
  dom_name_field = (config[:parsepfs][:acc]) ? 5 : 6
  consider_clans = config[:parsepfs][:clans]
  pid_regexp = config[:parsepfs][:PIDregexp]
  remove_empty = config[:parsepfs][:rempty]

  count = 0
  current_id = nil 
  last_id = nil 
  current_xdom = Array.new

  IO.foreach(result_file) {|line|
    next if (line[0] == '#')

    if m = pid_regexp.match(line)
      fields = line.chomp!.split
      current_id = m[1]
      if ( current_id != last_id and (not last_id.nil?) )
        count += 1
        print "\rParsing entry  %d... " % [count]
        if (config[:parsepfs][:resolve_overlaps])
          current_xdom = resolve_overlaps(current_xdom)
        end
        outfile.puts current_xdom.join("\n") unless (remove_empty and current_xdom.length == 1)
        last_id = current_id
        current_xdom = [">#{current_id}"]
      end 
      current_xdom << ">#{current_id}" if current_xdom.empty?
      name = fields[dom_name_field]
      if (consider_clans)
         name = fields[14] unless(fields[14] == 'No_clan')
      end
      next if ( (not config[:parsepfs][:cutoff].nil? ) and (fields[12].to_f > config[:parsepfs][:cutoff]) )
      current_xdom << "#{fields[dom_start_field]}\t#{fields[dom_end_field]}\t#{name}\t#{fields[12]}"
      last_id = current_id
    end 
  }
  if (config[:parsepfs][:resolve_overlaps])
    current_xdom = resolve_overlaps(current_xdom)
  end
  outfile.puts current_xdom.join("\n") unless (remove_empty and current_xdom.length == 1)
  outfile.close
  print "\rParsing entry  %d... done\n\n" % [count]
end

#### OVERLAP RESOLUTION ####
def resolve_overlaps(xdom_array)
  return xdom_array if (xdom_array.length == 2)
  prev_pos = 1
  current_pos = 2 
  loop do 
    break if (current_pos > xdom_array.length-1)
    current_line = xdom_array[current_pos].split
    prev_line = xdom_array[prev_pos].split
    if (current_line[0].to_i <= prev_line[1].to_i)
      (prev_line[3].to_f <= current_line[3].to_f) ? xdom_array.delete_at(current_pos) : xdom_array.delete_at(prev_pos)
      next
    end
    prev_pos += 1
    current_pos += 1
  end
  return xdom_array
end

### MAIN ###
def main
  if RUBY_VERSION.to_f < 1.9
    STDERR.puts "Requires >= Ruby 1.9 (running version #{RUBY_VERSION})"
    exit(-1)
  end
  config = read_config()
  outfiles = split_fasta(config)
  outfiles = run_pfamscan(outfiles, config)
  result_file = merge_results(outfiles, config)
  create_xdom(result_file, config) if (config[:files][:parse])
end

if ARGV.length == 1 then main() else puts "Usage: #{$0} <path-to-config.yaml>" end

