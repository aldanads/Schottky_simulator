# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 12:18:07 2022

@author: ALDANADS
"""

import cv2
import os

image_folder = 'Crystal growth'
video_name = 'Crystal.avi'

images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
# https://stackoverflow.com/questions/47294468/avoid-lexicographic-ordering-of-numerical-values-with-python-min-max
# Sort by numbers, not lexicographically
images = sorted(images, key = lambda x: int(x.split('_')[0]))
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape

video = cv2.VideoWriter(video_name, 0, 5, (width,height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()