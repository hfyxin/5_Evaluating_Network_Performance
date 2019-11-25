# Predict objects for all images in given folder, using yolo (Darknet) network. 
# Based on Ben's original notebook yoloV2analysisImportImages.ipynb
#
# Darknet package must be in the current directory:
# - darkflow (need to be built)
# - cfg (network config)
# - ckpt (network weights)
# Darknet tutorial: https://www.youtube.com/watch?v=fSM6cdFQdwI&t
#

from darkflow.net.build import TFNet
import cv2
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import sys

# I/O files
ImPath = '..\\All_Combined\\Validation\\'    # input image folder
filenames = []                               # image filenames
CSVpath = '..\\All_Combined\\Detected_CSV_Results\\'   # output results


# Configuration. See tutorial mentioned above
options = {
    'model': 'cfg/tiny-yolo-7c.cfg',
    'load': 29250,
    'threshold': 0.2,
    'labels': 'labelsAll.txt'
}
tfnet = TFNet(options)   # Display NN summary


# Get the list of filenames
for r, d, f in os.walk(ImPath): # r=root, d=directories, f = files
    for fn in f:
        if '.jpg' in fn:
            filenames.append(fn)
print("Total {:d} image files to be processed in folder: {}".format(len(filenames), ImPath))


# For displaying progress
def show_progress(count, total):
    # Display the progress in a fancy way.
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s/%s\r' % (bar, count, total))
    sys.stdout.flush()


# Warning
input("The output will erase existing files in folder: {}\n".format(CSVpath) +
    "Press any key to start, or Ctrl + C to stop:")
    
# Do the detection, save results to individual .csv files
i = 0
for fn in filenames:
    # Process image. See tutorial mentioned above
    img = cv2.imread(ImPath+fn, cv2.IMREAD_COLOR)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    results = tfnet.return_predict(img)
    
    # Save to csv
    with open(CSVpath+fn[:-4]+'.csv', 'w') as csvFile:
        fields = ['label','confidence','topleft','bottomright']
        writer = csv.DictWriter(csvFile,fieldnames=fields)      #,dialect='myDialect')
        writer.writeheader()
        writer.writerows(results)
        
    i += 1
    show_progress(i, len(filenames))

print('\n')