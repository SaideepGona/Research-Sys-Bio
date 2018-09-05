# This script is for emptying the cluster scratch drives programmatically

import os
import glob
import sys
import shutil

user_dir = sys.argv[1]
scratch_path = "/scratch/"
user_path = scratch_path + user_dir + "/"

print(os.listdir(user_path))
sys.stdout.flush()


shutil.rmtree(user_path)

