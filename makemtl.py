import sys
import os.path

if len(sys.argv) < 3:
  print "  Usage: " + sys.argv[0] +" model.obj addtomaterials.mtl"
  print "  Usage: for %i in (*.obj) do makemtl.py %i C:\Thesis\SDF\sdfDefault.mtl"
  sys.exit(2)

mtls = set()

def gather(file, keyword):
  for line in file:
    words = line.split();
    if len(words) > 0 and words[0] == keyword:
      mtls.add(words[1])


f = open(sys.argv[1], "r")
if (os.path.exists(sys.argv[2])):
  of = open(sys.argv[2], "r+")
  gather(of, "newmtl")
else:
  of = open(sys.argv[2], "w")

print len(mtls);
gather(f, "usemtl")
print len(mtls);

of.seek(0);
for colname in mtls:
  col = colname.split('_')
  r = float(col[1])/255
  g = float(col[2])/255
  b = float(col[3])/255
  of.write("newmtl " + colname + "\n")
  of.write("Ka " + str(r*0.15)[:7] + " " + str(g*0.15)[:7] + " " + str(b*0.15)[:7] + "\n")
  of.write("Kd " + str(r)[:7] + " " + str(g)[:7] + " " + str(b)[:7] + "\n")
  of.write("Ks " + str(r)[:7] + " " + str(g)[:7] + " " + str(b)[:7] + "\n")
  of.write("illum 2\nNs 23\n\n")
print("done.");