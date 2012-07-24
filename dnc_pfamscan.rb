#!/usr/bin/env ruby

WDIR        = '.'
FASTADIR    = "#{WDIR}/fas"
HMMOUTDIR   = "#{WDIR}/pfsouts"
MAXTHREADS  = 8


def split_fasta(opt, total)
  Dir.mkdir(FASTADIR) unless File.exists?(FASTADIR)
  fasta = File.new(opt, "r")
  outfiles = []
  fasta_basename = File.basename(opt, ".*") 
  in_file_no = 0
  out_file = File.new("#{FASTADIR}/#{fasta_basename}-#{outfiles.size}.fa", "w")
  IO.foreach(fasta) do |line|
    if (line[0] == '>')
      if (in_file_no == total)
        STDERR.puts "done."
        out_file.close()
        outfiles << out_file
        out_file = File.new("#{FASTADIR}/#{fasta_basename}-#{outfiles.size}.fa", "w")
        in_file_no = 0
      end
      in_file_no += 1
      STDERR.print "\rSplittig [part #{outfiles.size}] pfamseq [#{in_file_no}]... "
    end
    out_file.puts line
  end
  return outfiles
end


def run_pfamscan(files, bin)

  Dir.mkdir(HMMOUTDIR) unless File.exists?(HMMOUTDIR)


#  files.each do |fasta|
#    outfile_base = File.basename(fasta, ".fa")
#    puts "#{bin} -fasta #{fasta.path} -dir pfam_scan_db -clan_overlap -cpu 2 -outfile #{HMMOUTDIR}/#{outfile_base}.pfsout}"
#  end

end


def merge_results(opt)


end


def main

  pfam_scan_path  = ARGV[0]
  chunk_size      = ARGV[1]
  fasta_file      = ARGV[2]
  outfiles = split_fasta(pfam_scan_path, chunk_size.to_i)
  run_pfamscan(outfiles, pfam_scan_path)

end

if (ARGV.length == 3)
  main()
else
  STDERR.puts "Usage: #{$0} <path-to-pfamscan> <fasta-chunk-size> <fastaFile>"
end
