function Mr = registerfilm(M, FrameRef)

fig = figure(1);
imshow(M(:,:,FrameRef));
% [xmin ymin width height]
rect = floor(getrect(fig))

rectindX = [rect(1):(rect(1)+rect(3))];
rectindY = [rect(2):(rect(2)+rect(4))];

template = M(rectindY, rectindX,FrameRef);

bigrectX = [max(1,rect(1)-20):min(size(M,2),(rect(1)+rect(3)+20))];
bigrectY = [max(1,rect(2)-20):min(size(M,1),(rect(2)+rect(4)+20))];

% bigrectX=1:size(M,2);
% bigrectY=1:size(M,1);

Mr = M;

for i=1:size(M,3)
   
    out = normxcorr2(template,M(bigrectY,bigrectX,i));
%     imshow(out)
%     size(out)
%     pause
    
    [maxval,ind] = max(abs(out(:)))
    [Yp,Xp] = ind2sub(size(out),ind)
    % bottom right corner matters
    Yoffset = (Yp-size(template,1))-(rect(2)-bigrectY(1))
    Xoffset = (Xp-size(template,2))-(rect(1)-bigrectX(1))
    
    % [BW,xi,yi] = roipoly(M(:,:,1));

    se = translate(strel(1), [-Yoffset -Xoffset]);
    Mr(:,:,i) = imdilate(M(:,:,i),se);
end;