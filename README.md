# get_LC
Get light curves (LC) from K2 and TESS space missions.

This module is built to intake a pandas dataframe of targets with columns of 'tic' or 'epic' identifier, 'ra', and 'dec'. 

If you do not have identifiers for your targets, but rather just 'ra' and 'dec', you can use the catalog_queries.py to find your identifiers. 

Check out the catalog_queries_quickstart.ipynb tutorial for getting target IDs, and once you have your identifiers saved with the ra/dec coordinates, head over to the get_LC_quickstart.ipynb to download your LCs.
