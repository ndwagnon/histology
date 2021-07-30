%Name: Noah Wagnon
%Date: 04-21-2021
%Assignment: 582 Final Project
%Description: Script to perform detection and segmentation of
%lung cell nuclei in H&E stained images

clear variables
clc; close all;

%load testing image and take complemented grayscale
originalImage = imread('lungscc66.jpeg');
testImage = imcomplement(im2gray(originalImage));
[m,n] = size(testImage);
figure;imshow(originalImage); title('Input Image');
figure; imshow(testImage); title('Grayscaled and Complemented');

%apply top hat and bottom hat filter twice
se = strel('disk',30);
top = imtophat(testImage,se);
figure;imshow(top); title('Top Hat filtering');
bottom = imbothat(testImage,se);
figure; imshow(bottom); title('Bottom Hat filtering');
topBottomFiltered = testImage + (top - bottom);
figure; imshow(topBottomFiltered); title('1 Top/Bottom iteration');

%Second iteration
se = strel('disk',15);
top = imtophat(topBottomFiltered,se);
bottom = imbothat(topBottomFiltered,se);
figure; imshow(top-bottom); title('Tophat filter minus Bottomhat filter');
topBottomFiltered = topBottomFiltered + (top - bottom);
figure; imshow(topBottomFiltered); title('2 Top/Bottom iterations');

% calculate histogram of top/bottom filtered image
numPixels = numel(testImage);
[counts,n] = imhist(topBottomFiltered);
pixelCumSum = cumsum(counts);
foreground = 0.99 * numPixels;
background = 0.215 * numPixels;
low = find(pixelCumSum>background, 1, 'first');

%Extra Contrast Enhancement
enhancedImage = imadjust(topBottomFiltered, [low/255 1.0],[0.0 1.0],1.8);
figure; imshow(enhancedImage); title('Gamma corrected image');

%image binarization. Dynamically select a threshold value via Otsu's method
reshapedImage = reshape(enhancedImage, numPixels, 1);
reshapedImage = sort(reshapedImage);
thresh = graythresh(reshapedImage(numPixels * .5:end));
binarizedImage = imbinarize(enhancedImage,thresh);
figure; imshow(binarizedImage); title('Binarized Image');

%open the image to clear debris
se = strel('disk', 5);
openedImage = imopen(binarizedImage, se);
figure; imshow(openedImage); title('Opened binarized Image');

%Compute center points with connected components analysis
connected = bwconncomp(openedImage);
centers = regionprops(connected, 'Centroid');
size = connected.NumObjects;
centerMap = zeros(m,m);
for i = 1 : size
    centroid = round(centers(i).Centroid);
    centerMap(centroid(2), centroid(1)) = 1;
end

centerMap = centerMap > 0;

%dilate and display center map
figure; imshow(imdilate(centerMap, strel('disk', 3)));
title('Nuclei Center Map');
figure; imshow(originalImage + uint8(imdilate(centerMap * 255, strel('disk', 3))));
title('Center Map Overlaid on Original');

%Laplacian of Gaussian Edge filtering (size of 5, sigma of 3)
h = fspecial('log', 5, 3);
logFiltered = imfilter(double(openedImage)*255, h);
figure; imshow(logFiltered); title('After Log Filtering');

%Resize filtered output and add to original image
logFilteredResized = uint8(logFiltered .* (3^2));
finalImage = originalImage + (logFilteredResized*5);
figure; imshow(finalImage); title('Final Segmented Output');

%Display the final output with segmentation and detection
finalOutput = finalImage + uint8(imdilate(centerMap * 255, strel('disk',3)));
figure; imshow(finalOutput);
title('Final Output showing Segmentation and Center Points');