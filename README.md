# Robotics Publication Searching.

!CAUTION!
DBLP and Google Scholar does not shows every publications.
It only shows CS related publications.
- - -
For German Database
1. Write PI names and affiliations into 'PIlist.xlsx' 
2. Run 'pubsearch.py'
3. Copy the raw data 'publications.xlsx' to 'DBLPpub_post_processed.xlsx/database'.
4. Copy PI infos to 'DBLPpub_post_processed.xlsx/Affiliation', and sort Affiliation AtoZ (!IMPORTANT!).
5. Run 'duplication.py'
6. Copy the data of 'duplication.xlsx' to 'DBLPpub_post_processed.xlsx/Duplication check'.
7. Refresh 'DBLPpub_post_processed.xlsx/PI Summary' (Drag A5~U5 and copy-drag(?) the functions enough until every PIs listed + Set year range)
8. Check if 'DBLPpub_post_processed.xlsx/Institution Summary' shows every institutions (check the last institution).
- - -
For MIRMI Database
1. Write PI names and affiliations into 'PIlist_MIRMI.xlsx'
2. Run 'mirmi.py'
3. Copy the raw data 'publications.xlsx' to 'DBLPpub_post_processed_MIRMI.xlsx/database'.
4. Run 'scholar.py' (this is only for Nature & Science papers) *Sience Robotics is crawled by 'mirmi.py'
5. Copy the raw data 'publication_details.xlsx' to 'DBLPpub_post_processed_MIRMI.xlsx/database'.
6. Copy PI infos to 'DBLPpub_post_processed.xlsx/Affiliation'
7. Run 'MIRMI_jointpaper.py'
8. Copy the data of 'duplication.xlsx' to 'DBLPpub_post_processed_MIRMI.xlsx/Joint publications'.
9. Refresh 'DBLPpub_post_processed_MIRMI.xlsx/PI Summary' (Drag A5~O5 and copy-drag(?) the functions enough until every PIs listed + Set year range)
