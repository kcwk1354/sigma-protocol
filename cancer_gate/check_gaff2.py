import openmmforcefields, os, glob
pkg_dir = os.path.dirname(openmmforcefields.__file__)
gaff_xmls = glob.glob(pkg_dir + '/**/gaff-2.11.xml', recursive=True)
print('Found:', gaff_xmls)
if gaff_xmls:
    with open(gaff_xmls[0]) as f:
        lines = f.readlines()
    print('Total lines:', len(lines))
    targets = ['"c"', '"c3"', '"ca"', '"f"', '"n3"', '"o "', '"oh"', '"hc"', '"ha"', '"ho"']
    for line in lines[:300]:
        if any(t in line for t in targets):
            print(line.rstrip())
