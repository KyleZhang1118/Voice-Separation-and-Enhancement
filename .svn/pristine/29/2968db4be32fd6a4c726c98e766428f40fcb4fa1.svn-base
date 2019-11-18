function [Percent] = PrintLoopPCw(tt,NumFrames,oldPercent)
%PrintLoopPCw  Prints computation percentage on screen
%
% [PC] = PrintLoopPCw(LOOP,NUMLOOPS,oldPC) prints the current percentage of
% computations done in a loop. Use and initialize as in the following 
% examples. Also enables the use with WHILE loops.
%
%    PrintLoopPCw('\n\n   Starting computations. ');
%    for LOOP=1:NUMLOOPS,
%       PrintLoopPCw(LOOP,NUMLOOPS);
%
%       % ...
%       % ... your loop code ...
%       % ...
%    end
%
% Nested FOR loops example:
% 
%    PrintLoopPCw('\n\n   Starting computations. ');
%    for jj=1:NUMJJ,
%       for kk=1:NUMKK,
%          PrintLoopPCw((jj-1)*NUMKK+kk,NUMJJ*NUMKK);
%
%          % ... 
%          % ... your loop code ...
%          % ... 
%       end
%    end
%
% WHILE loop example:
%
%    LOOP = 1;
%    PC = PrintLoopPCw('\n\n   Starting computations. ');
%    while LOOP<=NUMLOOPS,
%       PC = PrintLoopPCw(LOOP,NUMLOOPS,PC);
%
%       % ... 
%       % ... your loop code including "LOOP = LOOP + 1;" at some stage (even conditionally) ...
%       % ... 
%    end

% Release date: August 2008
% Author: Eric A. Lehmann, Perth, Australia (www.eric-lehmann.com)
%
% Copyright (C) 2008 Eric A. Lehmann
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

Percent = 0;

if nargin==0,		% initialisation
    fprintf('Loop execution: done 0%%  ');

elseif nargin==1,	% customised initialisation
    if ~ischar(tt),
        error('Single input parameter must be a string');
    end
    fprintf([tt 'Loop execution: done 0%%  ']);
    
elseif nargin==2,   % use with FOR loop
    Percent = floor(100*tt/NumFrames);
    PercentOld = floor(100*(tt-1)/NumFrames);
    if PercentOld~=Percent,
        fprintf('\b\b\b\b%-4s',[num2str(Percent) '%']);
    end
    if tt==NumFrames,	% safe to say that this is the last call to PrintLoopPC...
        fprintf('\n');
    end

else        % use with WHILE loop
    Percent = floor(100*tt/NumFrames);
    if oldPercent~=Percent,
        fprintf('\b\b\b\b%-4s',[num2str(Percent) '%']);
        if tt==NumFrames,	% safe to say that this is the last call to PrintLoopPC...
            fprintf('\n');
        end
    end
end
