import xml.etree.ElementTree as ET
tree = ET.parse('./BlastUniref90.xml')
root = tree.getroot()

def main():
    
    arrayProteins=[]
    
# for hit in root:
#    print(hit,id)
#    print(description,description.attrib)
for hit in root.iter('hit'):
    print(hit.attrib)
    for querySeq in hit:
        print(querySeq.text)

if __name__ == '__main__':
    main()  