# Robotics Publication Searching.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dblpAPI is modified for the fast and reliable crawling.
Please use the library within this project. Not online installed dblpAPI library.

Install requirements
1. pandas
2. openpyxl
3. xlsxwriter
4. biopython

!CAUTION!
DBLP and Google Scholar does not shows every publications.
It only shows CS related publications.
- - -
For German Database
1. Write PI names and affiliations into 'PIlist.xlsx' !!Because of DBLP's strange searching mechanism, some names are not properly work. Please refer 'PIsearching key.xlsx' for the existing names. Please check new names in DBLP webpage.
3. Run 'pubsearch.py'
4. Copy the raw data 'publications.xlsx' to 'DBLPpub_post_processed.xlsx/database', and sort name AtoZ (!IMPORTANT!).
5. Copy PI infos to 'DBLPpub_post_processed.xlsx/Affiliation', and sort Affiliation AtoZ (!IMPORTANT!).
6. Run 'duplication.py' #IMPORTANT: Every Institutions that needed to be check should be in script! #IMPORTANT: Every Publication titles that needed to be check should be in script!
7. Copy the data of 'duplication.xlsx' to 'DBLPpub_post_processed.xlsx/Duplication check'. 
8. Refresh 'DBLPpub_post_processed.xlsx/PI Summary' (Drag A5~U5 and copy-drag(?) the functions enough until every PIs listed + Set year range)
9. Check if 'DBLPpub_post_processed.xlsx/Institution Summary' shows every institutions (check the last institution).
- - -
For MIRMI Database
1. Write PI names and affiliations into 'PIlist_MIRMI.xlsx' !!Because of DBLP's strange searching mechanism, some names are not properly work. Please refer 'PIsearching key.xlsx'
2. Run 'mirmi.py'
3. Copy the raw data 'publications.xlsx' to 'DBLPpub_post_processed_MIRMI.xlsx/database'.
4. Run 'scholar.py' (this is only for Nature & Science papers) *Sience Robotics is crawled by 'mirmi.py'
5. Copy the raw data 'publication_details.xlsx' to 'DBLPpub_post_processed_MIRMI.xlsx/database'.
6. Copy PI infos to 'DBLPpub_post_processed.xlsx/Affiliation'
7. Run 'MIRMI_jointpaper.py'
8. Copy the data of 'duplication.xlsx' to 'DBLPpub_post_processed_MIRMI.xlsx/Joint publications'.
9. Refresh 'DBLPpub_post_processed_MIRMI.xlsx/PI Summary' (Drag A5~O5 and copy-drag(?) the functions enough until every PIs listed + Set year range)
---
For pubmed: Not completed: duplication issues, journal title choosing issues
https://medium.com/@felipe.odorcyk/scrapping-data-from-pubmed-database-78a9b53de8ca

ToDo: Solve duplication problems

1. Run 'pubmed_test.py'.
2. Copy data from 'pubmed_result.xlsx' to 'DBLPpub_post_processed_pubmedadd.xlsx/database', and sort name AtoZ (!IMPORTANT!).
