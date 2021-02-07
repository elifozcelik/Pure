
# A function for getting the significant hit's (e-value) unique protein IDs in a list. 
def evalue_checker(text_name):
  # Getting the indexes of the hits part, lefting out the alignments. 
  idx1=0
  idx2=0
  file1=open(text_name)
  lines=file1.readlines()
  for i in range(0,len(lines)):
    if lines[i].find("Sequences producing significant alignments:") != -1:
      idx1=i
    elif lines[i].find(">") != -1:
      idx2=i
      break

  # idx1+2 & idx2-2 only gives the hits list
  only_hit_IDs=[]
  for i in range (idx1+2,idx2-2):
    one_list=lines[i].split(" ")
    filter_object = filter(lambda x: x != "", one_list)
    one_list1 = list(filter_object)
    # one_list1[0]==Uniprot ID, one_list last element-2 == e-value 
    length=len(one_list1)
    one_list1[length-1]=one_list1[length-1].strip()
    filter_object = filter(lambda x: x != "", one_list1)
    one_list2 = list(filter_object)
    length2=len(one_list2)
    print(one_list2)
    #print("**"+one_list1[length-2])
    if float(one_list2[length2-1])<1e-5: #using 0.00001
      only_hit_IDs.append(one_list2[0])

  return (only_hit_IDs)

def fastareader(filename):
  seqDict = {}
  header=""
  filein = open(filename, 'r')
  for line in filein:
    if line.startswith('>'):
      header = line[1:].strip()
      seqDict[header] = ""
    else:
      seqDict[header] += line.strip()
  return seqDict

def unique_file_creator(blast1,blast2,outfile):
  hit1=evalue_checker(blast1)
  hit2=evalue_checker(blast2)
  unique_IDS=[]

  for ID in hit1:
    if ID not in hit2:
      unique_IDS.append(ID)

  for ID in hit2:
    if ID not in unique_IDS:
      unique_IDS.append(ID)


  # Unique_IDs would have only unique hits, getting their sequences from database:

  db_dict=fastareader("all_eu.fasta")
  out = open(outfile, 'w') 

  for ID in unique_IDS:
    for name in db_dict.keys():
      if name.find(ID) != -1:
        out.write(">"+name+"\n")
        out.write(db_dict[name]+"\n")
        break
        
  out.close()
  return 0

# Putting the files that are going to be examined:
disease_blast_file_names=["disease_3_blast","disease7.out","disease9.out","disease12.out","disease13.out","disease15.out","disease20.out","disease22.out","disease24.out"]
ndisease_blast_file_names=["ndisease_3_blast","ndisease7.out","ndisease9.out","ndisease12.out","ndisease13.out","ndisease15.out","ndisease20.out","ndisease22.out","ndisease24.out"]
out_file_names=["unique_blast3","unique_blast7","unique_blast9","unique_blast12","unique_blast13","unique_blast15","unique_blast20","unique_blast22","unique_blast24"]
for i in range(0,9):
  unique_file_creator(disease_blast_file_names[i],ndisease_blast_file_names[i],out_file_names[i]) 
