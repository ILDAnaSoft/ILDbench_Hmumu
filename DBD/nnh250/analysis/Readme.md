# Caution
This part contains many DESY-specific part.
You need to be careful.

# How To Use
0. Assume you already initialized ilcsoft and include own library. If not, do it (written in ../mymarlin/Readme.md).
1. Prepare your list which contains the path of slcio file you want to proceed.
2. Edit `hmumu.xml` (check your parameters carefully).
3. Edit `job.sh` by hand, then do `. job.sh`. `tmp/` directory will be created, and it contains all xml files.
4. Edit `submit.sh` by hand, then do `. submit.sh`.
5. You should have processed root files at `output/PROCESS_NAME/SOMETHING.root`.

# Explanation
This place is used for job submission.


## Each File
- `list_*.txt`: external file which contain the path of slcio file you want to proceed
- `runMarlin`: prepare lots of xml files and job controlling files, this works when you do `. job.sh`
- `submit.sh`: for submitting all jobs
