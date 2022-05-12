function SCV=det_scv(DET)
    SCV=det_moment(DET,2)/det_moment(DET,1)^2-1;
end