# -*- coding: utf-8 -*-
"""database_seperation.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1CH5t_xho-AvO_utp6Z6dLgqzvA1GsuiG
"""

import json

gene_file= open("all_pages.txt",'r')
all=gene_file.read()
list1=all.split("\n")
group_dict={ }
key=""

# Creating the nested dictionary
for i in range(0,len(list1)): 
  if list1[i]=="":
    continue;
  if list1[i].find("Group") != -1:
    key_num=list1[i]
    group_dict[key_num]={}
  else:
    if list1[i].find("ENSG") != -1:
      key=list1[i]
      group_dict[key_num][key]=""
    else:
      group_dict[key_num][key]+=(list1[i]+"   ")

only_one_disease_gene_groups=[]
multiple_disease_gene_groups=[]
none_disease_gene_groups=[]
only_one_disease_gene_names=[]
    
for group in group_dict.keys():
  #equals 0 for each group
  count=0
  for gene in group_dict[group].keys():
    if(group_dict[group][gene].count("   ")>6):
      wanted_gene=gene
      count+=1
  
  if (count==1):
    only_one_disease_gene_groups.append(group)
    only_one_disease_gene_names.append(wanted_gene)
  if (count>1):
    multiple_disease_gene_groups.append(group)
  if(count==0):
    none_disease_gene_groups.append(group)

out = open('one_disease_groups.txt', 'w')   
for group in only_one_disease_gene_groups:
  out.write(group+"\n")
out.close()

multiple_different_disease_groups=0
multiple_different_disease_groups_names=[]

# For multiple disease groups, we have to check if the diseases are identical.
for group_name in multiple_disease_gene_groups:
  diseases=[]
  for gene_name in group_dict[group_name].keys():
    if(group_dict[group_name][gene_name].count("   ")>6):
      list2=group_dict[group_name][gene_name].split("   ")
      last=len(list2) #space
      for i in range(6,last): #must be the disease names
        diseases.append(list2[i])

  while(" " in diseases) : 
    diseases.remove(" ") 

  # Searching same diseases for each group:
  x=0
  for i in diseases:
    if(diseases.count(i)==1):
      x=1

  if (x):
    multiple_different_disease_groups+=1
    multiple_different_disease_groups_names.append(group_name)

out = open('multiple_different_gene_groups.txt', 'w')   
for group in multiple_different_disease_groups_names:
  out.write(group+"\n")
out.close()

#print("Groups:",len(group_dict.keys()))
#print("No disease groups:",len(none_disease_gene_groups))
#print("Only one disease groups:",len(only_one_disease_gene_groups))
#print("Multiple different disease groups:",multiple_different_disease_groups)

out = open('one_disease_genes_in_their_groups.txt', 'a')   
for name in only_one_disease_gene_names:
  out.write(name+"\n")
out.close()

id=1
groups=[]
for key in group_dict.keys():
  group_type=""
  if key in multiple_different_disease_groups_names:
    group_type="multiple disease genes"
  elif key in only_one_disease_gene_groups:
    group_type="only 1 disease gene"
  else:
    group_type="not interested"
  
  genes=[]
  for gene in group_dict[key]:
      list2=group_dict[key][gene].split("  ")
      while("" in list2) : 
        list2.remove("") 
      
      name=list2[0]
      description=list2[1]
      chromosome=list2[2]
      strand=list2[3]
      start_pos=list2[4]
      end_pos=list2[5]
      # If it has disease there will be 6 and more elements in list2

      gene_diseases=[]
      if len(list2)>6:
        for i in range(6,len(list2)):
          gene_diseases.append(list2[i])

      dict2={
          'gene_name':gene,
          'name':name,
          'description':description,
          'chromosome':chromosome,
          'strand':strand,
          'start_pos':start_pos,
          'end_pos':end_pos,
          'diseases':gene_diseases

      }
      genes.append(dict2)
      dict2={}

  dict1={'id':key,
         'type':group_type,
         'genes':genes
         }
  groups.append(dict1)
  dict1={}
  #id+=1


#Groups[0] == Group 1 | index
# i.values()= id'sini, tipini ve genlerinin olduğu listeyi içeriyor. i.keys() id,type,genes'e esit.
#ornek: multiple disease genesleri boyle cekebilirsin direk: 
#for i in groups:
  #if "only 1 disease gene" in i.values():
    #print(json.dumps(i, indent=2))

#!pip install pyensemblrest

from ensemblrest import EnsemblRest
from time import sleep
paralog_included_trees=[]
paralogs_not_in_tree=[]

# Displaying a specific genes group

for i in groups:
  if "only 1 disease gene" in i.values():
    gene_list= i["genes"]
    for a in gene_list:
      ensembl_id=a["gene_name"]   #a["gene_name"] --> ensembl gene ID'ye eşit.
      if ensembl_id== "ENSG00000116147":
        print(json.dumps(i, indent=2))


# Code for 30 trees that will examined normally 

for i in range (1,31):
  print("Disease gene",i,only_one_disease_gene_names[i])

# Getting trees from ensembl

for i in groups:
  if "only 1 disease gene" in i.values():
    gene_list= i["genes"]
    for a in gene_list:
      ensembl_id=a["gene_name"]
      if ensembl_id in only_one_disease_gene_names: #hastalık yapan genin agacina bakiyorum sadece

        if (ensembl_id=="ENSG00000203690"):
          continue

        tree=ensRest.getGeneTreeMemberById(id=ensembl_id,nh_format="simple", content_type="text/x-nh")
        human_count=tree.count("ENSP0") #insan sayısını sayıyorum
        if human_count>=2:
          paralog_included_trees.append(ensembl_id)
        else:
          paralogs_not_in_tree.append(ensembl_id)

out = open('paralog_included_trees_disease_genes.txt', 'a')   
for name in paralog_included_trees:
  out.write(name+"\n")
out.close()

out = open('paralogs_not_in_tree_disease_genes', 'a')   
for name in paralogs_not_in_tree:
  out.write(name+"\n")
out.close()

ensRest = EnsemblRest()
tree=ensRest.getGeneTreeMemberById(id='ENSG00000116985',nh_format="simple", content_type="text/x-nh")
print(json.dumps(groups[26], indent=2))
print(tree.count("ENSP0"))