function[msk] = morph_mask(I)

% Find regions surrounding 'double edges'

seB=strel('square',2); % Big 'diamond'????
seS=strel('square',2); % Small
% se1=strel('square',2);
msk=imdilate(I,seB);
msk=imerode(msk,seB);
msk=imerode(msk,seS);
msk=imdilate(msk,seS);
msk=logical(abs(msk-1));