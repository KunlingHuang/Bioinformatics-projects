Universe = vanilla
executable = /usr/bin/bash
Arguments = 2-1.hmagma.sh -c $(chr) -g $(gene) 
Request_memory = 1M
Request_cpus = 1
#concurrency_limits = khuang82:125
Log = ../condor_log/2-1/chr$(chr)/simple.$(gene).log
Output = ../condor_log/2-1/chr$(chr)/simple.$(gene).out
error = ../condor_log/2-1/chr$(chr)/simple.$(gene).err
# Send the job to Held state on failure.
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
periodic_release = (NumJobStarts < 10) && ((CurrentTime - EnteredCurrentStatus) > 300)

# Queue chr, gene, start, end, dir from ./condor_submit_txt/gene_chr.test.txt

# Queue chr, gene, start, end, dir from ./condor_submit_txt/gene_chr.location.chr1_chr10.txt

Queue chr, gene, start, end, dir from ./condor_submit_txt/gene_chr.location.chr11_chr22.txt


