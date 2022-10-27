import gcn_altered.scraper as scraper
#import shutil

scraperOBJ=scraper.Scraper() 
scraperOBJ.scrape()       #COLLECTS THE WHOLE GCN CIRCULAR ARCHIVE DATA AND UPDATES THE JSON FILE WITH ALL THE NEW UPDATES

#shutil.rmtree('/gcncc/data/gcn3') 
