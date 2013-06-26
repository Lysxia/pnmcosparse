function [ Ohm ] = FDAmatrix( d )
%FDAmatrix Finite Difference Analysis matrix

    Ohm = convmtx([-1 1],d-1);

end

