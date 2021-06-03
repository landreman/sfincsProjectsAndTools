#!/usr/bin/env python
import re

def filename2ref(fullfilename):
    s1 = fullfilename.rsplit('.',1)
    if not (s1[-1] == "bc"):
        raise Warning("File ending is not '.bc': " + fullfilename)
    filename = s1[0] # remove file ending if found
    s2 = filename.rsplit('/',1)
    if len(s2)>1:
        filename = filename.rsplit('/',1)[1] # remove dirs if applicable
        

    # try to see if equilibrium type is in filename
    try:
        i = filename.index("EIM")
        if i>=0:
            return "EIM"
    except ValueError:
        pass
    try:
        i = filename.index("KJM")
        if i>=0:
            return "KJM"
    except ValueError:
        pass
    
    # try to see if refXXX is in the filename
    ref = "ref"
    try:
        i = filename.index(ref)
        if i>=0:
            i = i + len(ref)
        return ref + re.findall(r'\d+', filename[i:])[0]
    except ValueError:
        pass

    # finally, just return the first string digits of suitable length (<4)
    digistrings = re.findall(r'\d+', filename)
    for ds in digistrings:
        if len(ds) < 4:
            return ref + ds

    # if we are here, we have failed
    raise ValueError("Could not deduce reference discharge number from filename for '" + fullfilename + "'.\nTry including 'ref' in front of the W7X reference number, or 'EIM' or 'KJM' directly in the filename.")

if __name__=="__main__":
    import sys
    argv = sys.argv
    argc = len(sys.argv)
    if argc>1:
        fn = argv[1]
        print(filename2ref(fn))
    else:
        print("Usage: ./filename2ref <.bc filename>")
