function [outdata,mind,maxd] = normalise(data, p2, p3)
indata = double(data);
percent = p3;
ndata = prod(size(indata));
[val,pos] = sort(indata(:));
upos = round(ndata*percent);
maxd = val(upos);
mind = min(indata(:));
outdata = (indata-double(mind)*ones(size(indata))) / (maxd-mind);

%outdata(find(outdata(:)>1)) = 1;
sz=size(outdata);
ind=find(outdata>1);

n1ind=ind+1;
n1ind(n1ind>sz(1)*sz(2)*sz(3))=sz(1)*sz(2)*sz(3);


n2ind=ind-1;
n2ind(n2ind<1)=1;


n3ind=ind+sz(1);

bname=floor(n3ind./(sz(1)*sz(2)));
bname2=floor(ind./(sz(1)*sz(2)));
miind=find(bname~=bname2);
n3ind(miind)=n1ind(miind);


n4ind=ind-sz(1);

bn=floor(n4ind./(sz(1)*sz(2)));
bn2=floor(ind./(sz(1)*sz(2)));
miind2=find(bn~=bn2);
n4ind(miind2)=n2ind(miind2);

outdata(outdata(n1ind)>1)=0.5;
outdata(outdata(n2ind)>1)=0.5;
outdata(outdata(n3ind)>1)=0.5;
outdata(outdata(n4ind)>1)=0.5;

outdata(ind)=(outdata(n1ind)+outdata(n2ind)+outdata(n3ind)+outdata(n4ind))./4;
