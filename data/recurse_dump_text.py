# recurse_dump_text.py - utility function to dump XML tags and text.
#
# Usage:
#   python recurse_dump_text.py file.xml
#
# Author: Randy Heiland
#

import xml.etree.ElementTree as ET
import sys
# import shutil

show_debug = False

argc = len(sys.argv)
if show_debug:
    print('argc=',argc)
if argc == 2:
    xml_file1 = sys.argv[1]
else:
    print('Usage: file1.xml')

# The original, hierarchical XML
if show_debug:
    print("\n\n     recurse_dump_text.py: ------- parsing ",xml_file1)
tree1 = ET.parse(xml_file1)
root1 = tree1.getroot()
for child in root1:
  print(child.tag)
  for child2 in child:
    print('--',child2.tag,child2.text)
    for child3 in child2:
      print('----',child3.tag,child3.text)
      for child4 in child3:
        print('------',child4.tag,child4.text)
#  print(child.tag, child.attrib)
