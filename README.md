# Robotics Publication Searching.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dblpAPI is modified for the fast and reliable crawling.
Please use the library within this project. Not online installed dblpAPI library.

Install requirements
1. pandas
2. openpyxl
3. xlsxwriter

!CAUTION!
DBLP and Google Scholar does not shows every publications.
It only shows CS related publications.
- - -
For German Database
1. Write PI names and affiliations into 'PIlist.xlsx' !!Because of DBLP's strange searching mechanism, some names are not properly work. Please refer 'PIsearching key.xlsx'
3. Run 'pubsearch.py'
4. Copy the raw data 'publications.xlsx' to 'DBLPpub_post_processed.xlsx/database'.
5. Copy PI infos to 'DBLPpub_post_processed.xlsx/Affiliation', and sort Affiliation AtoZ (!IMPORTANT!).
6. Run 'duplication.py'
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
For pubmed:
https://medium.com/@felipe.odorcyk/scrapping-data-from-pubmed-database-78a9b53de8ca
under construction
