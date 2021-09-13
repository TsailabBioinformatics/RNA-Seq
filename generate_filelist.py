"""
This program makes a file.list for the RNAseq pipeline
"""
import os

# list of file names
names = []
# list of file names after cleaning
names_cleaned = []

def store_paths(path):
    """
    Extracts and stores paths of fastq files into an array

    @param path: the base directory path
    """
    with os.scandir(path) as fastq:
        for element in fastq:
            if element.is_file() and element.name.endswith("fastq"):
                names.append(element.name)
            if element.is_dir():
                store_paths(element.path)

def clean_list():
    """
    Chops unimportant nomenclature 
    """
    for name in names:
        clean_name = ""
        i = 0
        while name[i] != '.':
            clean_name += name[i]
            i += 1
        names_cleaned.append(clean_name)

def write_list():
    """
    Writes all of the file names to a list (file.list)
    """
    with open("file.list", "wt") as list:
        for name in names_cleaned:
            list.write(name)
            list.write("\n")
        
store_paths("./fastq")
clean_list()
# removes duplicates
names_cleaned = list(set(names_cleaned))
write_list()
