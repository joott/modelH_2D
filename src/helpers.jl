import ParallelStencil: INDICES
ix, iy = INDICES[1], INDICES[2]

function NNp(n)
    n%L+1
end

function NNp2(n)
    (n+1)%L+1
end

function NNm(n)
    (n+L-2)%L+1
end

function NNm2(n)
    (n+L-3)%L+1
end

macro d_xc(A) esc(:( ($A[$NNp($ix),$iy] - $A[$NNm($ix),$iy])/2 )) end
macro d_yc(A) esc(:( ($A[$ix,$NNp($iy)] - $A[$ix,$NNm($iy)])/2 )) end

macro prd_d_xc(A,B) esc(:( ($A[$NNp($ix),$iy] * $B[$NNp($ix),$iy]
                          - $A[$NNm($ix),$iy] * $B[$NNm($ix),$iy])/2 )) end
macro prd_d_yc(A,B) esc(:( ($A[$ix,$NNp($iy)] * $B[$ix,$NNp($iy)]
                          - $A[$ix,$NNm($iy)] * $B[$ix,$NNm($iy)])/2 )) end

macro d2_xy(A)
    esc(:(
         ($A[$NNp2($ix),$iy] + $A[$NNm2($ix),$iy]
        + $A[$ix,$NNp2($iy)] + $A[$ix,$NNm2($iy)]
        - 4*$A[$ix,$iy]) * 0.25
    ))
end

function view_tuple(u)
    if size(u, 3) == 2
        return (@view(u[:,:,1]),@view(u[:,:,2]))
    end

    (@view(u[:,:,1]),@view(u[:,:,2]),@view(u[:,:,3]))
end
