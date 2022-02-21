# recurse_xml.py - utility function to convert a hierarchical PhysiCell config file (.xml), i.e.,
#                   with inheritance of parents in <cell_definitions>,
#                   to a "flattened" (expanded) config file.
#
# Usage:
#   python recurse_xml.py <hierarchical.xml> <flat.xml> 
#
# - args are optional; by default 1st= PhysiCell_settings.xml, 2nd= tmp2.xml
#
# Author: Randy Heiland
#

import xml.etree.ElementTree as ET
import sys
# import shutil

show_debug = False

# xml_hier_file = "PhysiCell_settings.xml"
# xml_flat_file = "tmp2.xml"
# xml_flat_file = "tmp_flat.xml"
argc = len(sys.argv)
if show_debug:
    print('argc=',argc)
if argc == 3:
    xml_file1 = sys.argv[1]
    xml_file2 = sys.argv[2]
else:
    print('Usage: file1.xml file2.xml')

# The original, hierarchical XML
if show_debug:
    print("\n\n     recurse_xml.py: ------- parsing ",xml_file1)
tree1 = ET.parse(xml_file1)
root1 = tree1.getroot()

tree2 = ET.parse(xml_file2)
root2 = tree2.getroot()

# # The newly constructed, expanded (flattened) XML
# # shutil.copy(xml_hier_file, xml_flat_file)
# tree_flat = ET.parse(xml_flat_file)
# root_flat = tree_flat.getroot()
# cell_defs_flat = root_flat.find("cell_definitions")

# num_nodes = 0
# #for node in root.iter():
# for node in []:
#   if show_debug:
#     print(node)
#   num_nodes += 1
# #print("num_nodes=",num_nodes)

# if show_debug:
#     print('----- recursive -----')
# prefix = ""
# full_path = ""
# last_tag = ""
# # root_cell_def = xml_root.find("cell_definitions//cell_definition")
# cell_defs_hier = root1.find("cell_definitions")

# #--------------------------
# def print_children(parent):
#     global cell_def_name, secretion_substrate_name, prefix, full_path, last_tag #, len_leaf_text
#     for child in parent:
#         if child.tag == 'cell_definition':
#             cell_def_name = child.attrib['name']
#             if show_debug:
#                 print("==== cell_def_name=",cell_def_name,'\n')
#         elif child.tag == 'substrate' and parent.tag =='secretion':  # i.e., //secretion//substrate
#         # elif child.tag == 'substrate':
#             if show_debug:
#                 print("\n===================== FOUND //secretion//substrate!! child.tag=",child.tag)
#                 print("\n===================== parent.tag=",parent.tag)
#                 print("===================== name=",child.attrib['name'])
#             secretion_substrate_name = child.attrib['name']
#             # print("===================== exiting!\n")
#             # sys.exit()

#         if show_debug:
#             print(prefix,child.tag,child.attrib)
#         full_path += '//' + child.tag
#         # full_path += '//' + child.tag
#         last_tag = child.tag
#         if show_debug:
#             print('full_path=',full_path)
#     # print(child.attrib)
#         prefix = prefix + "--"
#         print_children(child)
# #  prefix = prefix - "--"
#     #-----------------------------
#     # Reached a leaf element at this point. Update the corresponding element in the new, expanded XML.
#     # full_path = full_path[0:len(full_path)-len(last_tag)-2]
#     # print("   leaf path=",full_path[2:])
#     try:
#         leaf_text = cell_defs_hier.find(full_path[2:]).text
#         leaf_text = leaf_text.strip('\n')
#         leaf_text = leaf_text.strip('\t')
#         leaf_text = leaf_text.strip()
#         if show_debug:
#             print('    leaf_text=',leaf_text)
#         len_leaf_text = len(leaf_text)
#         if show_debug:
#             print('    len_leaf_text=',len_leaf_text)

#         # Update the expanded XML
#         leaf_text2 = cell_defs_flat.find(full_path[2:]).text
#         if show_debug:
#             print('    leaf_text2=',leaf_text2)
#         if (len_leaf_text > 0):
#             find_str = full_path[2:]
#             new_str = "cell_definition[@name='" + cell_def_name + "']"
#             if show_debug:
#                 print(">>>>>>>>>>>    new_str=",new_str)
#             find_str = find_str.replace("cell_definition",  new_str)
#             if find_str.find('//secretion//substrate//'):
#                 # <substrate name="director signal">
#                 if show_debug:
#                     print(">>>>>>>>>>>>>>>>>>   found //secretion//substrate//")
#                 new_substrate_str = "substrate[@name='" + secretion_substrate_name + "']"
#                 find_str = find_str.replace("substrate",  new_substrate_str)
#             if show_debug:
#                 print("  *** find_str=",find_str)
#                 print("  *** for cell_def_name ",cell_def_name,": replace ",cell_defs_flat.find(find_str).text," with ", leaf_text)
#             # cell_defs_flat.find(find_str).text = leaf_text
#             cell_defs_flat.find(find_str).text = cell_defs_hier.find(find_str).text
#     except:
#         pass
#     try:
#         idx = full_path.rindex("//")  # strip off the last element
#         full_path = full_path[0:idx]
#     except:
#         pass
#     if show_debug:
#         print("~~~ popped out: full_path=",full_path)
    
#     # root_cell_def = xml_root.find("cell_definitions + ")

#     prefix = prefix[2:]

# # print_children(root)
# print_children(cell_defs_hier)

# # new_xml_file = "tmp3.xml"
# # new_xml_file = "recurse_xml_out.xml"
# new_xml_file = "flat_xml_out.xml"
# print("---> ",new_xml_file)
# tree_flat.write(new_xml_file)

# sys.exit(1)

print('-----------------------')
for child in root1:
  print(child.tag)
  for child2 in child:
    print('--',child2.tag,child2.text)
    for child3 in child2:
      print('----',child3.tag,child3.text)
      for child4 in child3:
        print('------',child4.tag,child4.text)
#  print(child.tag, child.attrib)
