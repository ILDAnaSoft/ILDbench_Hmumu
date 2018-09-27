# Caution
This part contains many DESY-specific part.
You need to be careful.

# How To Use
0. Assume you already initialized ilcsoft and include own library. If not, do it (written in ../mymarlin/Readme.md).
1. Prepare your list which contains the path of slcio file you want to proceed.
2. Edit `hmumu.xml` (check your parameters carefully).
3. Do `. job.sh`.
4. Edit `submit.sh`, then do `. submit.sh`.
5. You should have processed root files as `output/PROCESS_NAME/SOMETHING.root`.

# Explanation
This place is used for job submission.


## Each File