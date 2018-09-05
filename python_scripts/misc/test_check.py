import sys
import subprocess

def QuickCheck(file):
    """
    Runs the samtools quickckeck command on a given input file to check for truncation. If there is a problem,
    will cancel the job and and write the filename to stdout
    """
    checkOut = subprocess.run(["samtools", "quickcheck", "-v", file])
    print(checkOut)

    if checkOut.returncode != 0:
        print("quickcheck failed")
        sys.stderr.write("Truncated File:" + file)
        sys.exit()

QuickCheck("truncated.bam")

print("not truncated")