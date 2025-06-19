function vertexassn = parcelassn_to_vertexassn(parcelassn,parceldata)
vertexassn = zeros(size(parceldata));
for iparcel = setdiff(unique(parceldata)',0)
    vertexassn(parceldata==iparcel) = parcelassn(iparcel);
end
end