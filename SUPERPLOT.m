function [] = SUPERPLOT(nRows, nCols, iRow, iCol, varargin)

%5th - 9th arguments
%sMar, tMar, bMar, cSpace, vSpace


if nargin<9
   vSpace = 0.02;
else
    vSpace = varargin{5};
end

if nargin<8
   cSpace = 0.02;
else
    cSpace = varargin{4};
end

if nargin<7
    sMar = .05;
else
    sMar = varargin{3};
end

 if nargin<6
    tMar = .05;
else
    tMar = varargin{2};
 end

 if nargin<5
    bMar = .1;
else
    bMar = varargin{1};
end

%% 

height = (1-tMar-bMar-(nRows-1)*vSpace)/(nRows);
width = (1-(nCols-1)*cSpace-2*sMar)/(nCols); 
vB = 1-(tMar+(iRow)*(height+vSpace));
hL = sMar+(iCol-1)*(width+cSpace);
subplot('Position',[hL(1), vB, width+(hL(end)-hL(1)), height]);
end
