from Bio import Entrez
import pandas as pd
import numpy as np
import re
from unidecode import unidecode

def search(query):
    Entrez.email = 'jaesug.jung@tum.de'
    handle = Entrez.esearch(db='pubmed', sort='relevance', retmax='250000', retmode='xml', term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'jaesug.jung@tum.de'
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
    results = Entrez.read(handle)
    return results

file_path = 'DBLPpub_post_processed.xlsx'
sheet_name = 'Affiliation'
data = pd.read_excel(file_path, sheet_name=sheet_name)
query = data.iloc[:, 0].tolist()
aff = data.iloc[:,1].tolist()

author_list=[]
title_list= []
abstract_list=[]
journal_list = []
language_list =[]
pubdate_year_list = []
pubdate_month_list = []
aff_list=[]

for k, q in enumerate(query):
    # q_tmp = unidecode(re.sub(r'[^a-zA-Z\s]', '', q))
    q_tmp = re.sub(r'[^a-zA-ZäöüÄÖÜß\s]', '', q)
    studies = search(q_tmp)

    # Check if the search result is not empty
    if 'IdList' not in studies or not studies['IdList']:
        print(f"Skipping {q_tmp} - No search results")
        continue
    
    studiesIdList = studies['IdList']

    studies = fetch_details(studiesIdList)
    chunk_size = 10000

    for chunk_i in range(0, len(studiesIdList), chunk_size):
        chunk = studiesIdList[chunk_i:chunk_i + chunk_size]
        papers = fetch_details(chunk)

        for i, paper in enumerate (papers['PubmedArticle']):
            author_list.append(q)
            title_list.append(paper['MedlineCitation']['Article']['ArticleTitle'])
            journal_list.append(paper['MedlineCitation']['Article']['Journal']['Title'])
            aff_list.append(aff[k])
            try:
                pubdate_year_list.append(int(paper['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']['Year']))
            except:
                pubdate_year_list.append('No Data')
    
    print(f'{q}, {aff[k]}, {i+1}')

df = pd.DataFrame(list(zip(author_list, pubdate_year_list, journal_list, title_list, aff_list)),columns=['Author','Year', 'Journal', 'Title','Affiliation'])

# Save DataFrame to Excel file
output_file_path = 'pubmed_result.xlsx'
df.to_excel(output_file_path, index=False)

print(f"Results saved to {output_file_path}")
