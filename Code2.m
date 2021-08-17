
%interval [0,180]

% angular_filter_bank(BW,fname)
%BW    - bandwidth (radians)
%fname - file name

%----------------------------------------------------------
function angular_filter_bank(BW,fname)
close all;
%---------------
%parameters
%---------------
FFTN    =   32;%64;%32;
TSTEPS  =   12; %15 degrees interval
DELTAT  =   pi/TSTEPS;
%---------------
%precompute
%---------------
[x,y]   =   meshgrid(-FFTN/2:FFTN/2-1,-FFTN/2:FFTN/2-1);
r       =   sqrt(x.^2+y.^2);
th      =   atan2(y,x);
th(th<0)=   th(th<0)+2*pi;  %unsigned
filter  =   [];

%-------------------------
%generate the filters
%-------------------------
for t0  =   0:DELTAT:(TSTEPS-1)*DELTAT
     t1     = t0+pi;                                %for the other lobe
     %-----------------
     %first lobe
     %-----------------
     d          = angular_distance(th,t0);
     msk        = 1+cos(d*pi/BW); 
     msk(d>BW)  = 0;
     rmsk       = msk;                              %save first lobe

     %-----------------
     %second lobe
     %-----------------
     d          = angular_distance(th,t1);
     msk        = 1+cos(d*pi/BW); 
     msk(d>BW)  = 0;
     rmsk       = (rmsk+msk);

     imagesc(rmsk);pause;
     rmsk   = transpose(rmsk);
     filter = [filter,rmsk(:)];
end;
     %-----------------
     %save the filters
     %-----------------
     angf  = filter;
     eval(sprintf('save %s angf',fname));

%----------------------
%write to a C file
%----------------------
fp = fopen(sprintf('%s.h',fname),'w');
fprintf(fp,'{\n');
for i = 1:size(filter,2)
    i
    k = 1;
    fprintf(fp,'{');
    for j = 1:size(filter,1)
        fprintf(fp,'%f,',filter(j,i));
        if(k == 32) k=0; fprintf(fp,'\n'); end;
        k = k+1;
    end;
    fprintf(fp,'},\n');
end;
fprintf(fp,'};\n');
fclose(fp);
%end function radial_filter_bank

%-----------------------------------------
%angular_distance
%computes angular distance-acute angle 
%-----------------------------------------
function d = angular_distance(th,t0)
    d = abs(th-t0);
    d = min(d,2*pi-d);
%end function angular_distance

%------------------------------------------------------------------------
%compute_coherence
%Computes the coherence image. 
%Usage:
%[cimg] = compute_coherence(oimg)
%oimg - orientation image
%cimg - coherence image(0-low coherence,1-high coherence)

function [cimg] = compute_coherence(oimg)
    [h,w]   =   size(oimg);
    cimg    =   zeros(h,w);
    N       =   2;
    %---------------
    %pad the image
    %---------------
    oimg    =   [flipud(oimg(1:N,:));oimg;flipud(oimg(h-N+1:h,:))]; %pad the rows
    oimg    =   [fliplr(oimg(:,1:N)),oimg,fliplr(oimg(:,w-N+1:w))]; %pad the cols
    %compute coherence
    for i=N+1:h+N
        for j = N+1:w+N
            th  = oimg(i,j);
            blk = oimg(i-N:i+N,j-N:j+N);
            cimg(i-N,j-N)=sum(sum(abs(cos(blk-th))))/((2*N+1).^2);
        end;
    end;
%end function compute_coherence

function n2 = dist2(x, c)


[ndata, dimx] = size(x);
[ncentres, dimc] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end

n2 = (ones(ncentres, 1) * sum((x.^2)', 1))' + ...
  		ones(ndata, 1) * sum((c.^2)',1) - ...
  		2.*(x*(c'));


function [enhimg, cimg, oimg,fimg,bwimg,eimg] =  fft_enhance_cubs(img, BLKSZ)
    global NFFT;

    if BLKSZ > 0
       NFFT        =   32;     %size of FFT
       OVRLP       =   2;      %size of overlap
       ALPHA       =   0.5;    %root filtering
       RMIN        =   4;%%3;      %min allowable ridge spacing
       RMAX        =   40;     %maximum allowable ridge spacing
       ESTRETCH    =   20;     %for contrast enhancement
       ETHRESH     =   19;      %threshold for the energy
    else
       NFFT        =   32;     %size of FFT
       BLKSZ       =   12;     %size of the block
       OVRLP       =   6;      %size of overlap
       ALPHA       =   0.5;    %root filtering
       RMIN        =   3;      %min allowable ridge spacing
       RMAX        =   18;     %maximum allowable ridge spacing
       ESTRETCH    =   20;     %for contrast enhancement
       ETHRESH     =   6;      %threshold for the energy
    end
    
    [nHt,nWt]   =   size(img);  
    img         =   double(img);    %convert to DOUBLE
    nBlkHt      =   floor((nHt-2*OVRLP)/BLKSZ);
    nBlkWt      =   floor((nWt-2*OVRLP)/BLKSZ);
    fftSrc      =   zeros(nBlkHt*nBlkWt,NFFT*NFFT); %stores FFT
    nWndSz      =   BLKSZ+2*OVRLP; %size of analysis window. 
    %-------------------------
    %allocate outputs
    %-------------------------
    oimg        =   zeros(nBlkHt,nBlkWt);
    fimg        =   zeros(nBlkHt,nBlkWt);
    bwimg       =   zeros(nBlkHt,nBlkWt);
    eimg        =   zeros(nBlkHt,nBlkWt);
    enhimg      =   zeros(nHt,nWt);
    
    %-------------------------
    %precomputations
    %-------------------------
    [x,y]       =   meshgrid(0:nWndSz-1,0:nWndSz-1);
    dMult       =   (-1).^(x+y); %used to center the FFT
    [x,y]       =   meshgrid(-NFFT/2:NFFT/2-1,-NFFT/2:NFFT/2-1);
    r           =   sqrt(x.^2+y.^2)+eps;
    th          =   atan2(y,x);
    th(th<0)    =   th(th<0)+pi;
    w           =   raised_cosine_window(BLKSZ,OVRLP); %spectral window

    %-------------------------
    %Load filters
    %-------------------------
    load angular_filters_pi_4;   %now angf_pi_4 has filter coefficients
    angf_pi_4 = angf;
    load angular_filters_pi_2;   %now angf_pi_2 has filter coefficients
    angf_pi_2 = angf;
    %-------------------------
    %Bandpass filter
    %-------------------------
    FLOW        =   NFFT/RMAX;
    FHIGH       =   NFFT/RMIN;
    
    dRLow       =   1./(1+(r/FHIGH).^4);    %low pass butterworth filter
    dRHigh      =   1./(1+(FLOW./r).^4);    %high pass butterworth filter
    dBPass      =   dRLow.*dRHigh;          %bandpass
    
    %-------------------------
    %FFT Analysis
    %-------------------------
    for i = 0:nBlkHt-1
        nRow = i*BLKSZ+OVRLP+1;  
        for j = 0:nBlkWt-1
            nCol = j*BLKSZ+OVRLP+1;
            %extract local block
            blk     =   img(nRow-OVRLP:nRow+BLKSZ+OVRLP-1,nCol-OVRLP:nCol+BLKSZ+OVRLP-1);
            %remove dc
            dAvg    =   sum(sum(blk))/(nWndSz*nWndSz);
            blk     =   blk-dAvg;   %remove DC content
            blk     =   blk.*w;     %multiply by spectral window
            %--------------------------
            %do pre filtering
            %--------------------------
            blkfft  =   fft2(blk.*dMult,NFFT,NFFT);
            blkfft  =   blkfft.*dBPass;             %band pass filtering
            dEnergy =   abs(blkfft).^2;
            blkfft  =   blkfft.*sqrt(dEnergy);      %root filtering(for diffusion)
            fftSrc(nBlkWt*i+j+1,:) = transpose(blkfft(:));
            dEnergy =   abs(blkfft).^2;             %----REDUCE THIS COMPUTATION----
            %--------------------------
            %compute statistics
            %--------------------------
            dTotal          =   sum(sum(dEnergy))/(NFFT*NFFT);
            fimg(i+1,j+1)   =   NFFT/(compute_mean_frequency(dEnergy,r)+eps); %ridge separation
            oimg(i+1,j+1)   =   compute_mean_angle(dEnergy,th);         %ridge angle
            eimg(i+1,j+1)   =   log(dTotal+eps);                        %used for segmentation
        end;%for j
    end;%for i

    %-------------------------
    %precomputations
    %-------------------------
    [x,y]       =   meshgrid(-NFFT/2:NFFT/2-1,-NFFT/2:NFFT/2-1);
    dMult       =   (-1).^(x+y); %used to center the FFT

    %-------------------------
    %process the resulting maps
    %-------------------------
    for i = 1:3
        oimg = smoothen_orientation_image(oimg);            %smoothen orientation image
    end;
    fimg    =   smoothen_frequency_image(fimg,RMIN,RMAX,5); %diffuse frequency image
    cimg    =   compute_coherence(oimg);                    %coherence image for bandwidth
    bwimg   =   get_angular_bw_image(cimg);                 %QUANTIZED bandwidth image
    %-------------------------
    %FFT reconstruction
    %-------------------------
    for i = 0:nBlkHt-1
        for j = 0:nBlkWt-1
            nRow = i*BLKSZ+OVRLP+1;            
            nCol = j*BLKSZ+OVRLP+1;
            %--------------------------
            %apply the filters
            %--------------------------
            blkfft  =   reshape(transpose(fftSrc(nBlkWt*i+j+1,:)),NFFT,NFFT);
            %--------------------------
            %reconstruction
            %--------------------------
            af      =   get_angular_filter(oimg(i+1,j+1),bwimg(i+1,j+1),angf_pi_4,angf_pi_2);
            blkfft  =   blkfft.*(af); 
            blk     =   real(ifft2(blkfft).*dMult);
            enhimg(nRow:nRow+BLKSZ-1,nCol:nCol+BLKSZ-1)=blk(OVRLP+1:OVRLP+BLKSZ,OVRLP+1:OVRLP+BLKSZ);
        end;%for j
    end;%for i
    %end block processing
    %--------------------------
    %contrast enhancement
    %--------------------------
    enhimg =sqrt(abs(enhimg)).*sign(enhimg);
    mx     =max(max(enhimg));
    mn     =min(min(enhimg));
    enhimg =uint8((enhimg-mn)/(mx-mn)*254+1);
    
    %--------------------------
    %clean up the image
    %--------------------------
    emsk  = imresize(eimg,[nHt,nWt]);
    enhimg(emsk<ETHRESH) = 128;
%end function fft_enhance_cubs

%-----------------------------------
%raised_cosine
%returns 1D spectral window
%syntax:
%y = raised_cosine(nBlkSz,nOvrlp)
%y      - [OUT] 1D raised cosine function
%nBlkSz - [IN]  the window is constant here
%nOvrlp - [IN]  the window has transition here
%-----------------------------------
function y = raised_cosine(nBlkSz,nOvrlp)
    nWndSz  =   (nBlkSz+2*nOvrlp);
    x       =   abs(-nWndSz/2:nWndSz/2-1);
    y       =   0.5*(cos(pi*(x-nBlkSz/2)/nOvrlp)+1);
    y(abs(x)<nBlkSz/2)=1;
%end function raised_cosine

%-----------------------------------
%raised_cosine_window
%returns 2D spectral window
%syntax:
%w = raised_cosine_window(blksz,ovrlp)
%w      - [OUT] 1D raised cosine function
%nBlkSz - [IN]  the window is constant here
%nOvrlp - [IN]  the window has transition here
%-----------------------------------
function w = raised_cosine_window(blksz,ovrlp)
    y = raised_cosine(blksz,ovrlp);
    w = y(:)*y(:)';
%end function raised_cosine_window

%---------------------------------------------------------------------
%get_angular_filter
%generates an angular filter centered around 'th' and with bandwidth 'bw'
%the filters in angf_xx are precomputed using angular_filter_bank.m
%syntax:
%r = get_angular_filter(t0,bw)
%r - [OUT] angular filter of size NFFTxNFFT
%t0- mean angle (obtained from orientation image)
%bw- angular bandwidth(obtained from bandwidth image)
%angf_xx - precomputed filters (using angular_filter_bank.m)
%-----------------------------------------------------------------------
function r = get_angular_filter(t0,bw,angf_pi_4,angf_pi_2)
    global NFFT;
    TSTEPS = size(angf_pi_4,2);
    DELTAT = pi/TSTEPS;
    %get the closest filter
    i      = floor((t0+DELTAT/2)/DELTAT);
    i      = mod(i,TSTEPS)+1; 
    if(bw == pi/4)
        r      = reshape(angf_pi_4(:,i),NFFT,NFFT)';
    elseif(bw == pi/2)
        r      = reshape(angf_pi_2(:,i),NFFT,NFFT)';
    else
        r      = ones(NFFT,NFFT);
    end;
%end function get_angular_filter


%-----------------------------------------------------------
%get_angular_bw_image
%the bandwidth allocation is currently based on heuristics
%(domain knowledge :)). 
%syntax:
%bwimg = get_angular_bw_image(c)
%-----------------------------------------------------------
function bwimg = get_angular_bw_image(c)
    bwimg   =   zeros(size(c));
    bwimg(:,:)    = pi/2;                       %med bw
    bwimg(c<=0.7) = pi;                         %high bw
    bwimg(c>=0.9) = pi/4;                       %low bw
%end function get_angular_bw


%-----------------------------------------------------------
%get_angular_bw_image
%the bandwidth allocation is currently based on heuristics
%(domain knowledge :)). 
%syntax:
%bwimg = get_angular_bw_image(c)
%-----------------------------------------------------------
function mth = compute_mean_angle(dEnergy,th)
    global NFFT;
    sth         =   sin(2*th);
    cth         =   cos(2*th);
    num         =   sum(sum(dEnergy.*sth));
    den         =   sum(sum(dEnergy.*cth));
    mth         =   0.5*atan2(num,den);
    if(mth <0)
        mth = mth+pi;
    end;
%end function compute_mean_angle

%-----------------------------------------------------------
%get_angular_bw_image
%the bandwidth allocation is currently based on heuristics
%(domain knowledge :)). 
%syntax:
%bwimg = get_angular_bw_image(c)
%-----------------------------------------------------------
function mr = compute_mean_frequency(dEnergy,r)
    global NFFT;
    num         =   sum(sum(dEnergy.*r));
    den         =   sum(sum(dEnergy));
    mr          =   num/(den+eps);
%end function compute_mean_angle

% FREQEST - Estimate fingerprint ridge frequency within image block
%
% Function to estimate the fingerprint ridge frequency within a small block
% of a fingerprint image.  This function is used by RIDGEFREQ
%
% Usage:
%  freqim =  freqest(im, orientim, windsze, minWaveLength, maxWaveLength)
%
% Arguments:
%         im       - Image block to be processed.
%         orientim - Ridge orientation image of image block.
%         windsze  - Window length used to identify peaks. This should be
%                    an odd integer, say 3 or 5.
%         minWaveLength,  maxWaveLength - Minimum and maximum ridge
%                     wavelengths, in pixels, considered acceptable.
% 
% Returns:
%         freqim    - An image block the same size as im with all values
%                     set to the estimated ridge spatial frequency.  If a
%                     ridge frequency cannot be found, or cannot be found
%                     within the limits set by min and max Wavlength
%                     freqim is set to zeros.
%
% Suggested parameters for a 500dpi fingerprint image
%   freqim = freqest(im,orientim, 5, 5, 15);
%
% See also:  RIDGEFREQ, RIDGEORIENT, RIDGESEGMENT
%
% Note I am not entirely satisfied with the output of this function.

% Peter Kovesi 
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% January 2005

    
function freqim =  freqest(im, orientim, windsze, minWaveLength, maxWaveLength)
    
    debug = 0;
    
    [rows,cols] = size(im);
    
    % Find mean orientation within the block. This is done by averaging the
    % sines and cosines of the doubled angles before reconstructing the
    % angle again.  This avoids wraparound problems at the origin.
    orientim = 2*orientim(:);    
    cosorient = mean(cos(orientim));
    sinorient = mean(sin(orientim));    
    orient = atan2(sinorient,cosorient)/2;

    % Rotate the image block so that the ridges are vertical
    rotim = imrotate(im,orient/pi*180+90,'nearest', 'crop');
    
    % Now crop the image so that the rotated image does not contain any
    % invalid regions.  This prevents the projection down the columns
    % from being mucked up.
    cropsze = fix(rows/sqrt(2)); offset = fix((rows-cropsze)/2);
    rotim = rotim(offset:offset+cropsze, offset:offset+cropsze);

    % Sum down the columns to get a projection of the grey values down
    % the ridges.
    proj = sum(rotim);
    
    % Find peaks in projected grey values by performing a greyscale
    % dilation and then finding where the dilation equals the original
    % values. 
    dilation = ordfilt2(proj, windsze, ones(1,windsze));
    maxpts = (dilation == proj) & (proj > mean(proj));
    maxind = find(maxpts);

    % Determine the spatial frequency of the ridges by divinding the
    % distance between the 1st and last peaks by the (No of peaks-1). If no
    % peaks are detected, or the wavelength is outside the allowed bounds,
    % the frequency image is set to 0
    if length(maxind) < 2
	freqim = zeros(size(im));
    else
	NoOfPeaks = length(maxind);
	waveLength = (maxind(end)-maxind(1))/(NoOfPeaks-1);
	if waveLength > minWaveLength & waveLength < maxWaveLength
	    freqim = 1/waveLength * ones(size(im));
	else
	    freqim = zeros(size(im));
	end
    end

    
    if debug
	%show(im,1)
	%show(rotim,2);
	figure(3),    plot(proj), hold on
	meanproj = mean(proj)
	if length(maxind) < 2
	    fprintf('No peaks found\n');
	else
	    plot(maxind,dilation(maxind),'r*'), hold off
	    waveLength = (maxind(end)-maxind(1))/(NoOfPeaks-1);
	end
    end
    
% NORMALISE - Normalises image values to 0-1, or to desired mean and variance
%
% Usage:
%             n = normalise(im)
%
% Offsets and rescales image so that the minimum value is 0
% and the maximum value is 1.  Result is returned in n.  If the image is
% colour the image is converted to HSV and the value/intensity component
% is normalised to 0-1 before being converted back to RGB.
%
%
%             n = normalise(im, reqmean, reqvar)
%
% Arguments:  im      - A grey-level input image.
%             reqmean - The required mean value of the image.
%             reqvar  - The required variance of the image.
%
% Offsets and rescales image so that it has mean reqmean and variance
% reqvar.  Colour images cannot be normalised in this manner.

% Copyright (c) 1996-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% January 2005 - modified to allow desired mean and variance


function n = normalise(im, reqmean, reqvar)

    if ~(nargin == 1 | nargin == 3)
       error('No of arguments must be 1 or 3');
    end
    
    if nargin == 1   % Normalise 0 - 1
	if ndims(im) == 3         % Assume colour image 
	    hsv = rgb2hsv(im);
	    v = hsv(:,:,3);
	    v = v - min(v(:));    % Just normalise value component
	    v = v/max(v(:));
	    hsv(:,:,3) = v;
	    n = hsv2rgb(hsv);
	else                      % Assume greyscale 
	    if ~isa(im,'double'), im = double(im); end
	    n = im - min(im(:));
	    n = n/max(n(:));
	end
	
    else  % Normalise to desired mean and variance
	
	if ndims(im) == 3         % colour image?
	    error('cannot normalise colour image to desired mean and variance');
	end

	if ~isa(im,'double'), im = double(im); end	
	im = im - mean(im(:));    
	im = im/std(im(:));      % Zero mean, unit std dev

	n = reqmean + im*sqrt(reqvar);
    end
    
%-----Sub functions-------
function [j,X, Y] = p(img, x, y, i)
% get pixel value based on chart:
%  4 | 3 | 2
%  5 |   | 1
%  6 | 7 | 8

switch (i)
    case {1, 9}
        Y=y;
        X=x+1;
        j = img(y, x + 1);
    case 2
        Y=y-1;
        X=x+1;
        j = img(y - 1, x + 1);
    case 3
        Y=y-1;
        X=x;
        j = img(y - 1, x);
    case 4
        Y=y-1;
        X=x-1;
        j = img(y - 1, x - 1);
    case 5
        Y=y;
        X=x-1;
        j = img(y, x - 1);
    case 6
        Y=y+1;
        X=x-1; 
        j = img(y + 1, x - 1);
    case 7
        Y=y+1;
        X=x;
        j = img(y + 1, x);
    case 8
        Y=y+1;
        X=x+1; 
        j = img(y + 1, x + 1);
end

% RIDGEFILTER - enhances fingerprint image via oriented filters
%
% Function to enhance fingerprint image via oriented filters
%
% Usage:
%  newim =  ridgefilter(im, orientim, freqim, kx, ky, showfilter)
%
% Arguments:
%         im       - Image to be processed.
%         orientim - Ridge orientation image, obtained from RIDGEORIENT.
%         freqim   - Ridge frequency image, obtained from RIDGEFREQ.
%         kx, ky   - Scale factors specifying the filter sigma relative
%                    to the wavelength of the filter.  This is done so
%                    that the shapes of the filters are invariant to the
%                    scale.  kx controls the sigma in the x direction
%                    which is along the filter, and hence controls the
%                    bandwidth of the filter.  ky controls the sigma
%                    across the filter and hence controls the
%                    orientational selectivity of the filter. A value of
%                    0.5 for both kx and ky is a good starting point.
%         showfilter - An optional flag 0/1.  When set an image of the
%                      largest scale filter is displayed for inspection.
% 
% Returns:
%         newim    - The enhanced image
%
% See also: RIDGEORIENT, RIDGEFREQ, RIDGESEGMENT

% Reference: 
% Hong, L., Wan, Y., and Jain, A. K. Fingerprint image enhancement:
% Algorithm and performance evaluation. IEEE Transactions on Pattern
% Analysis and Machine Intelligence 20, 8 (1998), 777 789.

% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% January 2005

function newim = ridgefilter(im, orient, freq, kx, ky, showfilter)

    if nargin == 5
        showfilter = 0;
    end
    
    angleInc = 3;  % Fixed angle increment between filter orientations in
                   % degrees. This should divide evenly into 180
    
    im = double(im);
    [rows, cols] = size(im);
    newim = zeros(rows,cols);
    
    [validr,validc] = find(freq > 0);  % find where there is valid frequency data.
    ind = sub2ind([rows,cols], validr, validc);

    % Round the array of frequencies to the nearest 0.01 to reduce the
    % number of distinct frequencies we have to deal with.
    freq(ind) = round(freq(ind)*100)/100;
    
    % Generate an array of the distinct frequencies present in the array
    % freq 
    unfreq = unique(freq(ind)); 
    
    % Generate a table, given the frequency value multiplied by 100 to obtain
    % an integer index, returns the index within the unfreq array that it
    % corresponds to
    freqindex = ones(100,1);
    for k = 1:length(unfreq)
        freqindex(round(unfreq(k)*100)) = k;
    end
    
    % Generate filters corresponding to these distinct frequencies and
    % orientations in 'angleInc' increments.
    filter = cell(length(unfreq),180/angleInc);
    sze = zeros(length(unfreq),1);
    
    for k = 1:length(unfreq)
        sigmax = 1/unfreq(k)*kx;
        sigmay = 1/unfreq(k)*ky;
        
        sze(k) = round(3*max(sigmax,sigmay));
        [x,y] = meshgrid(-sze(k):sze(k));
        reffilter = exp(-(x.^2/sigmax^2 + y.^2/sigmay^2)/2)...
                .*cos(2*pi*unfreq(k)*x);

        % Generate rotated versions of the filter.  Note orientation
        % image provides orientation *along* the ridges, hence +90
        % degrees, and imrotate requires angles +ve anticlockwise, hence
        % the minus sign.
        for o = 1:180/angleInc
            filter{k,o} = imrotate(reffilter,-(o*angleInc+90),'bilinear','crop'); 
        end
    end

%    if showfilter % Display largest scale filter for inspection
%        figure(7), imshow(filter{1,end},[]); title('filter'); 
%    end
    
    % Find indices of matrix points greater than maxsze from the image
    % boundary
    maxsze = sze(1);    
    finalind = find(validr>maxsze & validr<rows-maxsze & ...
                    validc>maxsze & validc<cols-maxsze);
    
    % Convert orientation matrix values from radians to an index value
    % that corresponds to round(degrees/angleInc)
    maxorientindex = round(180/angleInc);
    orientindex = round(orient/pi*180/angleInc);
    i = find(orientindex < 1);   orientindex(i) = orientindex(i)+maxorientindex;
    i = find(orientindex > maxorientindex); 
    orientindex(i) = orientindex(i)-maxorientindex; 

    % Finally do the filtering
    for k = 1:length(finalind)
        r = validr(finalind(k));
        c = validc(finalind(k));

        % find filter corresponding to freq(r,c)
        filterindex = freqindex(round(freq(r,c)*100));
        
        s = sze(filterindex);   
        newim(r,c) = sum(sum(im(r-s:r+s, c-s:c+s).*filter{filterindex,orientindex(r,c)}));
    end

    
% RIDGEFREQ - Calculates a ridge frequency image
%
% Function to estimate the fingerprint ridge frequency across a
% fingerprint image. This is done by considering blocks of the image and
% determining a ridgecount within each block by a call to FREQEST.
%
% Usage:
%  [freqim, medianfreq] =  ridgefreq(im, mask, orientim, blksze, windsze, ...
%                                    minWaveLength, maxWaveLength)
%
% Arguments:
%         im       - Image to be processed.
%         mask     - Mask defining ridge regions (obtained from RIDGESEGMENT)
%         orientim - Ridge orientation image (obtained from RIDGORIENT)
%         blksze   - Size of image block to use (say 32) 
%         windsze  - Window length used to identify peaks. This should be
%                    an odd integer, say 3 or 5.
%         minWaveLength,  maxWaveLength - Minimum and maximum ridge
%                     wavelengths, in pixels, considered acceptable.
% 
% Returns:
%         freqim     - An image  the same size as im with  values set to
%                      the estimated ridge spatial frequency within each
%                      image block.  If a  ridge frequency cannot be
%                      found within a block, or cannot be found within the
%                      limits set by min and max Wavlength freqim is set
%                      to zeros within that block.
%         medianfreq - Median frequency value evaluated over all the
%                      valid regions of the image.
%
% Suggested parameters for a 500dpi fingerprint image
%   [freqim, medianfreq] = ridgefreq(im,orientim, 32, 5, 5, 15);
%
% I seem to find that the median frequency value is more useful as an
% input to RIDGEFILTER than the more detailed freqim.  This is possibly
% due to deficiencies in FREQEST.
%
% See also: RIDGEORIENT, FREQEST, RIDGESEGMENT

% Reference: 
% Hong, L., Wan, Y., and Jain, A. K. Fingerprint image enhancement:
% Algorithm and performance evaluation. IEEE Transactions on Pattern
% Analysis and Machine Intelligence 20, 8 (1998), 777 789.

% Peter Kovesi  
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% January 2005

function [freq, medianfreq] = ridgefreq(im, mask, orient, blksze, windsze, ...
                                        minWaveLength, maxWaveLength)     
    
    [rows, cols] = size(im);
    freq = zeros(size(im));
    
    for r = 1:blksze:rows-blksze
        for c = 1:blksze:cols-blksze
          blkim = im(r:r+blksze-1, c:c+blksze-1);   
          blkor = orient(r:r+blksze-1, c:c+blksze-1);       
          
          freq(r:r+blksze-1,c:c+blksze-1) =  ...
              freqest(blkim, blkor, windsze, minWaveLength, maxWaveLength);
        end
    end

    % Mask out frequencies calculated for non ridge regions
    freq = freq.*mask;
    
    % Find median freqency over all the valid regions of the image.
    medianfreq = median(freq(find(freq>0)));  

% RIDGEORIENT - Estimates the local orientation of ridges in a fingerprint
%
% Usage:  [orientim, reliability] = ridgeorientation(im, gradientsigma,...
%                                             blocksigma, ...
%                                             orientsmoothsigma)
%
% Arguments:  im                - A normalised input image.
%             gradientsigma     - Sigma of the derivative of Gaussian
%                                 used to compute image gradients.
%             blocksigma        - Sigma of the Gaussian weighting used to
%                                 sum the gradient moments.
%             orientsmoothsigma - Sigma of the Gaussian used to smooth
%                                 the final orientation vector field.
% 
% Returns:    orientim          - The orientation image in radians.
%                                 Orientation values are +ve clockwise
%                                 and give the direction *along* the
%                                 ridges.
%             reliability       - Measure of the reliability of the
%                                 orientation measure.  This is a value
%                                 between 0 and 1. I think a value above
%                                 about 0.5 can be considered 'reliable'.
%
%
% With a fingerprint image at a 'standard' resolution of 500dpi suggested
% parameter values might be:
%
%    [orientim, reliability] = ridgeorient(im, 1, 3, 3);
%
% See also: RIDGESEGMENT, RIDGEFREQ, RIDGEFILTER

% Original version by Raymond Thai,  May 2003
% Reworked by Peter Kovesi           January 2005
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk


function [orientim, reliability] = ...
             ridgeorient(im, gradientsigma, blocksigma, orientsmoothsigma)
        
    [rows,cols] = size(im);
    
    % Calculate image gradients.
    sze = fix(6*gradientsigma);   if ~mod(sze,2); sze = sze+1; end
    f = fspecial('gaussian', sze, gradientsigma); % Generate Gaussian filter.
    [fx,fy] = gradient(f);                        % Gradient of Gausian.
    
    Gx = filter2(fx, im); % Gradient of the image in x
    Gy = filter2(fy, im); % ... and y
    
    % Estimate the local ridge orientation at each point by finding the
    % principal axis of variation in the image gradients.
   
    Gxx = Gx.^2;       % Covariance data for the image gradients
    Gxy = Gx.*Gy;
    Gyy = Gy.^2;
    
    % Now smooth the covariance data to perform a weighted summation of the
    % data.
    sze = fix(6*blocksigma);   if ~mod(sze,2); sze = sze+1; end    
    f = fspecial('gaussian', sze, blocksigma);
    Gxx = filter2(f, Gxx); 
    Gxy = 2*filter2(f, Gxy);
    Gyy = filter2(f, Gyy);
    
    % Analytic solution of principal direction
    denom = sqrt(Gxy.^2 + (Gxx - Gyy).^2) + eps;
    sin2theta = Gxy./denom;            % Sine and cosine of doubled angles
    cos2theta = (Gxx-Gyy)./denom;
       
    sze = fix(6*orientsmoothsigma);   if ~mod(sze,2); sze = sze+1; end    
    f = fspecial('gaussian', sze, orientsmoothsigma);    
    cos2theta = filter2(f, cos2theta); % Smoothed sine and cosine of
    sin2theta = filter2(f, sin2theta); % doubled angles
    
    orientim = pi/2 + atan2(sin2theta,cos2theta)/2;

    % Calculate 'reliability' of orientation data.  Here we calculate the
    % area moment of inertia about the orientation axis found (this will
    % be the minimum inertia) and an axis  perpendicular (which will be
    % the maximum inertia).  The reliability measure is given by
    % 1.0-min_inertia/max_inertia.  The reasoning being that if the ratio
    % of the minimum to maximum inertia is close to one we have little
    % orientation information. 
    
    Imin = (Gyy+Gxx)/2 - (Gxx-Gyy).*cos2theta/2 - Gxy.*sin2theta/2;
    Imax = Gyy+Gxx - Imin;
    
    reliability = 1 - Imin./(Imax+.001);
    
    % Finally mask reliability to exclude regions where the denominator
    % in the orientation calculation above was small.  Here I have set
    % the value to 0.001, adjust this if you feel the need
    reliability = reliability.*(denom>.001);


% RIDGESEGMENT - Normalises fingerprint image and segments ridge region
%
% Function identifies ridge regions of a fingerprint image and returns a
% mask identifying this region.  It also normalises the intesity values of
% the image so that the ridge regions have zero mean, unit standard
% deviation.
%
% This function breaks the image up into blocks of size blksze x blksze and
% evaluates the standard deviation in each region.  If the standard
% deviation is above the threshold it is deemed part of the fingerprint.
% Note that the image is normalised to have zero mean, unit standard
% deviation prior to performing this process so that the threshold you
% specify is relative to a unit standard deviation.
%
% Usage:   [normim, mask, maskind] = ridgesegment(im, blksze, thresh)
%
% Arguments:   im     - Fingerprint image to be segmented.
%              blksze - Block size over which the the standard
%                       deviation is determined (try a value of 16).
%              thresh - Threshold of standard deviation to decide if a
%                       block is a ridge region (Try a value 0.1 - 0.2)
%
% Returns:     normim - Image where the ridge regions are renormalised to
%                       have zero mean, unit standard deviation.
%              mask   - Mask indicating ridge-like regions of the image, 
%                       0 for non ridge regions, 1 for ridge regions.
%              maskind - Vector of indices of locations within the mask. 
%
% Suggested values for a 500dpi fingerprint image:
%
%   [normim, mask, maskind] = ridgesegment(im, 16, 0.1)
%
% See also: RIDGEORIENT, RIDGEFREQ, RIDGEFILTER

% Peter Kovesi         
% School of Computer Science & Software Engineering
% The University of Western Australia
% pk at csse uwa edu au
% http://www.csse.uwa.edu.au/~pk
%
% January 2005


function [normim, mask, maskind] = ridgesegment(im, blksze, thresh)
    
    im = normalise(im,0,1);  % normalise to have zero mean, unit std dev
    
    fun = inline('std(x(:))*ones(size(x))');
    
    stddevim = blkproc(im, [blksze blksze], fun);
    
    mask = stddevim > thresh;
    maskind = find(mask);
    
    % Renormalise image so that the *ridge regions* have zero mean, unit
    % standard deviation.
    im = im - mean(im(maskind));
    normim = im/std(im(maskind));    


%------------------------------------------------------------------------
%smoothen_frequency_image
%smoothens the frequency image through a process of diffusion
%Usage:
%new_oimg = smoothen_frequency_image(fimg,RLOW,RHIGH,diff_cycles)
%fimg       - frequency image image
%nimg       - filtered frequency image
%RLOW       - lowest allowed ridge separation
%RHIGH      - highest allowed ridge separation
%diff_cyles - number of diffusion cycles
%Contact:
%   ssc5@eng.buffalo.edu
%   www.eng.buffalo.edu/~ssc5
%Reference:
%S. Chikkerur, C.Wu and V. Govindaraju, "Systematic approach for feature
%extraction in Fingerprint Images", ICBA 2004
%------------------------------------------------------------------------
function nfimg = smoothen_frequency_image(fimg,RLOW,RHIGH,diff_cycles)
    valid_nbrs  =   3; %uses only pixels with more then valid_nbrs for diffusion
    [ht,wt]     =   size(fimg);
    nfimg       =   fimg;
    N           =   1;
    
    %---------------------------------
    %perform diffusion
    %---------------------------------
    h           =   fspecial('gaussian',2*N+1);
    cycles      =   0;
    invalid_cnt = sum(sum(fimg<RLOW | fimg>RHIGH));
    while((invalid_cnt>0 &cycles < diff_cycles) | cycles < diff_cycles)
        %---------------
        %pad the image
        %---------------
        fimg    =   [flipud(fimg(1:N,:));fimg;flipud(fimg(ht-N+1:ht,:))]; %pad the rows
        fimg    =   [fliplr(fimg(:,1:N)),fimg,fliplr(fimg(:,wt-N+1:wt))]; %pad the cols
        %---------------
        %perform diffusion
        %---------------
        for i=N+1:ht+N
         for j = N+1:wt+N
                blk = fimg(i-N:i+N,j-N:j+N);
                msk = (blk>=RLOW & blk<=RHIGH);
                if(sum(sum(msk))>=valid_nbrs)
                    blk           =blk.*msk;
                    nfimg(i-N,j-N)=sum(sum(blk.*h))/sum(sum(h.*msk));
                else
                    nfimg(i-N,j-N)=-1; %invalid value
                end;
         end;
        end;
        %---------------
        %prepare for next iteration
        %---------------
        fimg        =   nfimg;
        invalid_cnt =   sum(sum(fimg<RLOW | fimg>RHIGH));
        cycles      =   cycles+1;
    end;
%     cycles
%end function smoothen_orientation_image


%------------------------------------------------------------------------
%smoothen_orientation_image
%smoothens the orientation image through vectorial gaussian filtering
%Usage:
%new_oimg = smoothen_orientation_image(oimg)
%oimg     - orientation image
%new_oimg - filtered orientation image
%Contact:
%   ssc5@eng.buffalo.edu
%   www.eng.buffalo.edu/~ssc5
%Reference:
%M. Kaas and A. Witkin, "Analyzing oriented patterns", Computer Vision
%Graphics Image Processing 37(4), pp 362--385, 1987
%------------------------------------------------------------------------
function noimg = smoothen_orientation_image(oimg)
    %---------------------------
    %smoothen the image
    %---------------------------
    gx      =   cos(2*oimg);
    gy      =   sin(2*oimg);
    
    msk     =   fspecial('gaussian',5);
    gfx     =   imfilter(gx,msk,'symmetric','same');
    gfy     =   imfilter(gy,msk,'symmetric','same');
    noimg   =   atan2(gfy,gfx);
    noimg(noimg<0) = noimg(noimg<0)+2*pi;
    noimg   =   0.5*noimg;
%end function smoothen_orientation_image



% Returns if the single ridge in a bifurcation is 
% the going the same direction as the orientation (res=3)
% or in in the opposite (res=2)

function [res, progress, sx, sy, angle] = test_bifurcation(img, x, y, o, core_x, core_y)

   iax = 0; iay = 0; ibx = 0; iby = 0; icx = 0; icy = 0;
   progress = 1;
   path_len = 4; 

   pax = 0; pay=0; pbx = 0; pby=0; pcx = 0; pcy=0;
   pao = 0; pbo = 0; pco = 0;      %1-8 position of 3x3 square around minutia
   
   for i = 1:8
         [ta, xa, ya] = p(img, x, y, i);
         [tb, xb, yb] = p(img, x, y, i+1);
         if (ta > tb) 
           if pao == 0
              if i < 5
                 pao = 4 + i;
              else
                 pao = mod(4 + i, 9) + 1;
              end
              pax = xa;
              pay = ya;
           else
              if pbo == 0 
                 if i < 5
                    pbo = 4 + i;
                 else
                    pbo = mod(4 + i, 9) + 1;
                 end
                 pbx = xa;
                 pby = ya;
              else
                 if i < 5 
                    pco = 4 + i;
                 else
                    pco = mod(4 + i, 9) + 1;
                 end
                 pcx = xa;
                 pcy = ya;
                 break
              end
           end   
         end
   end

   xaa=0; yaa=0; xbb=0; ybb=0; xcc=0; ycc=0;
   
   hist=[];
   hist(1,:)=[pay, pax];
   hist(2,:)=[pby, pbx];
   hist(3,:)=[pcy, pcx];
   hist(4,:)=[y, x];
   stop=0;
   while ( progress < path_len && ~stop)
      progress = progress + 1;
      da = 0; db = 0; dc = 0;
      if pao ~= 0
          i = 1;

          cn = 0;
          for ii = 1:8
               [t1, x_A, y_A] = p(img, pax, pay, ii);
               [t2, x_B, y_B] = p(img,pax,pay, ii+1);
               cn = cn + abs (t1-t2);
          end
          cn = cn / 2;
 
          if cn == 1 || cn == 3
             stop=1;
          end

          while  i < 9 && da == 0

             [ta, xa, ya] = p(img, pax, pay, i);
             [tz, xz, yz] = p(img, pax, pay, i+1);

             ind_y=find(hist(:,1) == ya);
             ind_x=find(hist(ind_y,2) == xa);
             if numel(ind_x) > 0
                i = i+1;
                continue
             end

             if ta > tz && (xa ~=x || xa~=y)
                pax = xa;
                pay = ya;
                hist(size(hist,1)+1,:)=[pay, pax];
                da = 1;
                xaa = xa;
                yaa = ya;
             end
             i = i+1;
 
          end
          if da == 0 
             break
          end
      end

      if pbo ~= 0 && ~stop

          cn = 0;
          for ii = 1:8
               [t1, x_A, y_A] = p(img,pbx, pby, ii);
               [t2, x_B, y_B] = p(img,pbx, pby, ii+1);
               cn = cn + abs (t1-t2);
          end
          cn = cn / 2;
 
          if cn == 1 || cn == 3
             stop=1;
          end

          i=1;

          while  i < 9 && db == 0
 
             [ta, xa, ya] = p(img, pbx, pby, i);
             [tz, xz, yz] = p(img, pbx, pby, i+1);

             ind_y=find(hist(:,1) == ya);
             ind_x=find(hist(ind_y,2) == xa);
             if numel(ind_x) > 0
                i = i+1;
                continue
             end

             if ta > tz && (xa ~=x || xa~=y)
                pbx = xa;
                pby = ya;
                hist(size(hist,1)+1,:)=[pby, pbx];
                db=1;
                xbb=xa;                                                            
                ybb=ya;                                                            
             end
             i = i+1;
 
          end
      end

      if pco ~= 0 && ~stop

          cn = 0;
          for ii = 1:8
               [t1, x_A, y_A] = p(img,pcx, pcy, ii);
               [t2, x_B, y_B] = p(img, pcx, pcy, ii+1);
               cn = cn + abs (t1-t2);
          end
          cn = cn / 2;
 
          if cn == 1 || cn == 3
             stop=1;
          end


          i = 1;
          while  i < 9 && dc == 0
  
             [ta, xa, ya] = p(img, pcx, pcy, i);
             [tz, xz, yz] = p(img, pcx, pcy, i+1);

             ind_y=find(hist(:,1) == ya);
             ind_x=find(hist(ind_y,2) == xa);
             if numel(ind_x) > 0
                i = i+1;
                continue
             end

             if ta > tz && (xa ~=x || xa~=y)
                pcx = xa;
                pcy = ya;
                hist(size(hist,1)+1,:)=[pcy, pcx];
                dc = 1;
                xcc = xa;                                                            
                ycc = ya;                                                            
             end
             i = i+1;
           end
      end
   end

   d1 = sqrt(dist2([xaa yaa], [xbb ybb]));
   d2 = sqrt(dist2([xaa yaa], [xcc ycc]));
   d3 = sqrt(dist2([xcc ycc], [xbb ybb]));


   if d1 >= d3 &&  d2 >= d3
     sx = xaa;
     sy = yaa;
     ind = pao;
   elseif d1 >= d2 && d3 >= d2
     sx = xbb;
     sy = ybb;
     ind = pbo;
   elseif d3 >= d1 && d2 >= d1   
     sx = xcc;
     sy = ycc;
     ind = pco;
   else
     pause
   end


   d1 = sqrt(dist2([xaa yaa], [core_x core_y]));
   d2 = sqrt(dist2([xbb ybb], [core_x core_y]));
   d3 = sqrt(dist2([xcc ycc], [core_x core_y]));

   qx=0;
   qy=0;
   diff=0;
   if d1 >= d2 &&  d1 >= d3
     qx = xaa;
     qy = yaa;
     ind = pao;
   elseif d2 >= d3 && d2 >= d1
     qx = xbb;
     qy = ybb;
     ind = pbo;
   elseif d3 >= d2 && d3 >= d1
     qx = xcc;
     qy = ycc;
     ind = pco;
   else
     pause
   end


%sy
%sx
%y
%x
%   if core_x == 0
      angle = mod(atan2(y-sy, sx-x), 2*pi);
%      o
      if qx==sx && qy==sy %abs(angle - o) > pi/3
         res = 3;
      else 
         res = 3;
      end
%progress
%pause

%   else
%      if dist2([sy sx], [core_y core_x]) <= dist2([y x], [core_y core_x])
%         res = 2;
%      else 
%         res=3;
%      end 
%   end
















