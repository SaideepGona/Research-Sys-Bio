# This script is designed to be run alongside jobs running on the cluster. It periodically checks the size of
# stderr and stdout output files for unusually large file sizes as a stop-gap measure against unregulated behavior.
import os
import subprocess
import sys
import time

size_cut_off = 1000000000                    # File size which should not be exceeded

def size_check():
    all_files = os.listdir()
    std_files = []
    folder_size = 0

    for file in all_files:
        size = os.path.getsize(file)
        folder_size += size
        try:
            if file.endswith(".err") or file.endswith(".out"):
                if size > size_cut_off:
                    cancel_jobs()
                    sys.exit()
        except:
            continue

    return folder_size

def cancel_jobs():       # Cancels all jobs
    subprocess.run(["scancel", "-u", "sgona"])

while True:                           # While loop with sleep interval between SizeCheck executions
    total_sum = size_check()
    time.sleep(30)
