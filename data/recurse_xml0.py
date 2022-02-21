import xml.etree.ElementTree as ET
import sys

#tree = ET.parse('nanobio_settings.xml')
#tree = ET.parse('test1.xml')
tree = ET.parse('/Users/heiland/git/mcds2isa/HUVEC/HUVEC_v4_SHF_test.xml')
root = tree.getroot()

num_nodes = 0
#for node in root.iter():
for node in []:
  print(node)
  num_nodes += 1
#print("num_nodes=",num_nodes)

print('----- recursive -----')
prefix = ""
def print_children(parent):
  global prefix
  for child in parent:
    print(prefix,child.tag)
    prefix = prefix + "--"
    print_children(child)
#  prefix = prefix - "--"
  prefix = prefix[2:]

print_children(root)

sys.exit(1)

print('-----------------------')
for child in root:
  print(child.tag)
  for child2 in child:
    print('--',child2.tag)
    for child3 in child2:
      print('----',child3.tag)
      for child4 in child3:
        print('------',child4.tag)
#  print(child.tag, child.attrib)
