
%%
% Today I will show you how I plotted this treble clef in MATLAB:
%
% <<https://blogs.mathworks.com/steve/files/treble-clef-screen-shot.png>>
%
% My discovery and implementation process for doing this involved Unicode
% characters, SVG files, string processing, Bezier curves, |kron|, implicit
% expansion, and |polyshape|. I thought it was fun, and I learned a few
% things, so I wanted to share it.
%
% As a personal project, I make specialized fingering charts for the French
% horn. I do this in MATLAB, of course. For one of these charts, I wanted
% to be able to place several treble clef symbols, precisely and with high
% quality. That's I got started on this quest.
%
% First, I wondered if there is a Unicode symbol for treble clef (also
% known as the G clef). There is. One online resource I like to use for
% Unicode is <https://www.fileformat.info/info/unicode/index.htm
% FileFormat.info>, which has this info page for the treble clef:
%
% <<https://blogs.mathworks.com/steve/files/fileformat-info-screen-shot.png>>
%
% But I wasn't sure that using the Unicode character directly was really
% going to work out for me. The precise positioning and size of a Unicode
% character is too dependent on which font one uses. But the line
% <https://www.fileformat.info/info/unicode/char/1d11e/musical_symbol_g_clef.svg
% "Outline (as SVG file)"> on that info page caught my eye. Maybe I could
% use that?
% 
% I don't know that much about SVG files, but I downloaded the file and
% looked at it in a text editor. There didn't seem to be much to it:
%
% <<https://blogs.mathworks.com/steve/files/clef-svg-file-screen-shot.png>>
%
% Clearly, I needed to understand more about an SVG _path_. I found this
% <https://www.w3.org/TR/SVG/paths.html w3.org> page that explains SVG
% paths in detail. Fortunately for me, since I just wanted to try a quick
% prototype, the path in the treble clef outline file contained only a few
% commands:
%
% * "M" - move to a position
% * "L" - draw a line to a position
% * "Q" - draw a quadratic Bezier curve segment using three points
% * "Z" - close the current curve
%
% I certainly didn't want to write a fully general SVG path parser just to
% do a prototype. Instead, I used some MATLAB string processing
% functionality to break the path string down into easily processable
% chunks. Here's the string I was dealing with.
cubic=1;
load wrench_svg
clef_path=wrench{1};
% https://www.fileformat.info/info/unicode/char/1d11e/musical_symbol_g_clef.svg

%%
% I started with putting each command on a separate line.
commands = 'CMQLZ';
for k = 1:length(commands)
    ck = string(commands(k));
    clef_path = replace(clef_path,ck,newline + ck);
end

%%
% And then I split up the string by line breaks to make a string array.

clef_path = splitlines(clef_path);
% clef_path(1:15)

%%
% Now I need to explain about those "Q" commands. They indicate the drawing
% of a quadratic Bezier curve using three control points (the current point
% and the two additional points listed following the "Q" character). I
% wasn't familiar with the mathematics of drawing Bezier curves, but I
% quickly found <https://blogs.mathworks.com/graphics/2014/10/13/bezier-curves/ 
% Mike Garrity's old graphics blog post> on the subject. Here's a brief
% snippet of explanation and code from that old post.

P1 = [ 5; -10];
P2 = [18;  18];
P3 = [45;  15];

cla
placelabel(P1,'P_1');
placelabel(P2,'P_2');
placelabel(P3,'P_3');
xlim([0 50])
axis equal

%%
% _[Mike's explanation]_ And we've still got $t$ going from 0 to 1, but
% we'll use second-order polynomials like this:
%
% $$P(t) = (1-t)^2 P_1 + 2 (1-t) t P_2 + t^2 P_3$$
%
% These are known as the <http://en.wikipedia.org/wiki/Bernstein_polynomial
% Bernstein basis polynomials>.
%
% We can use the kron function to evaluate them.
t = linspace(0,1,101);
P = kron((1-t).^2,P1) + kron(2*(1-t).*t,P2) + kron(t.^2,P3);

hold on
plot(P(1,:),P(2,:))
hold off

%%
% Notice that the resulting curve starts at $P_1$ and ends at $P_3$. In
% between, it moves towards, but doesn't reach, $P_2$. _[End of Mike's
% explanation]_
%
% The |kron| function, used by Mike, computes something called the
% Kronecker tensor product. However, since implicit expansion was
% introduced in MATLAB three years ago, this is more easily computed simply
% by using element-wise multiplication.

t = t;
P = ((1-t).^2 .* P1) + (2*(1-t) .* t .* P2) + (t.^2 .* P3);

%%
% So, that's the computation I need to perform for each "Q" command in the
% SVG path. I don't think I need 101 points to draw each Bezier curve
% segment, though; let's go with 21. Also, I'm going to change the
% orientation to be a column.

t = linspace(0,1,21)';

%%
% Now let's build up |x| and |y| vectors based on each SVG path command.
if cubic
    l=4;
elseif quadratic
    l=3;
else
    l=1;
end
x = zeros(0,l);
y = zeros(0,l);

x_current = 0;
y_current = 0;
curves={};
ccounter=1;
cccounter=1;
currentcurve=zeros(0,1);
xxcounter=1;
for k = 1:length(clef_path)
    if clef_path(k) == ""
        continue
    end
    
    % Get the command character and a vector of the numeric values after
    % the command.
    command = extractBefore(clef_path(k),2);
    remainder = extractAfter(clef_path(k),1);
    remainder = strrep(remainder, ',', ' ')
    values = str2double(split(remainder));
    if isnan(values(1))
        values=values(2:(end-1));
    end
   
    switch command
        case "M"
            % Move to a position. Put a NaN-valued point in the vector.
            curves{ccounter}=currentcurve;
            currentcurve=zeros(0,l);
            cccounter=1;
            ccounter=ccounter+1;
            
            x = [x ; NaN];
            y = [y ; NaN];
            x_current = values(1);
            y_current = values(2);
            x(end+1,1) = x_current;
            y(end+1,1) = y_current;
%             curve{k}=[x;y];
        case "L"
            % Draw a line segment from the current point to the specified
            % point.
            if l==3
                currentcurve(cccounter:(cccounter+2),:)=[[x_current (x_current+values(1))/2 values(1)]; [y_current (y_current+values(2))/2 values(2)]; [1 1 1]];
            elseif l==4
                currentcurve(cccounter:(cccounter+2),:)=[[x_current (2/3*x_current+1/3*values(1)) (1/3*x_current+2/3*values(1)) values(1)]; [y_current (2/3*y_current+1/3*values(2)) (1/3*y_current+2/3*values(2)) values(2)]; [1 1 1 1]];
            end
            cccounter=cccounter+3;
            x_current = values(1);
            y_current = values(2);
            x(end+1,1) = x_current;
            y(end+1,1) = y_current;
%             curve{k}=[values x;y];
        case "Q"
            % Draw a quadratic Bezier curve segment using the current point
            % and the two additional points as control points.
            
            pt1 = [x_current y_current];
            pt2 = [values(1) values(2)];
            pt3 = [values(3) values(4)];
            currentcurve(cccounter:(cccounter+2),:)=[[x_current values(1) values(3)]; [y_current values(2) values(4)]; [1 1 1]];
            cccounter=cccounter+3;
            pts = ((1-t).^2 .* pt1) + (2*(1-t).*t .* pt2) + ...
                (t.^2 * pt3);
            x = [x ; pts(:,1)];
            y = [y ; pts(:,2)];
            
            x_current = values(3);
            y_current = values(4);
        case "C"
            % Draw a quadratic Bezier curve segment using the current point
            % and the two additional points as control points.
            
            pt1 = [x_current y_current];
            pt2 = [values(1) values(2)];
            pt3 = [values(3) values(4)];
            pt4 = [values(5) values(6)];
            currentcurve(cccounter:(cccounter+2),:)=[[x_current values(1) values(3) values(5)]; [y_current values(2) values(4) values(6)]; [1 1 1 1]];
            cccounter=cccounter+3;
            pts = ((1-t).^3 .* pt1) + (3*(1-t).^2.*t .* pt2) + ...
                (3*(1-t).*t.^2 * pt3) + ((t.^3.*pt4));
            x = [x ; pts(:,1)];
            y = [y ; pts(:,2)];
            
            x_current = values(5);
            y_current = values(6);
        case "Z"
            % TODO: curve-closing logic
            curves{ccounter}=currentcurve;
    end
    plot(x,y);
%     close
    end
curves{ccounter}=currentcurve;
%%
% Now, let's take a big breath and see what we've got!

plot(x,y)

%%
% *So close!* But upside down and oddly squished. Both those problems are
% easily fixed.

axis equal
axis ij

%%
% *Even closer!* But how do we get a filled shape?? That's where
% |polyshape| comes in. We've got |x| and |y| vectors containing four
% segments separated by |NaN| values. The four segments are for the outer
% curve, plus the outlines of the three holes or voids in the treble clef.
% The constructor for |polyshape| knows exactly how to deal with that kind
% of data. It can figure out automatically which curves are outer
% boundaries and which curves are inner boundaries. The resulting
% |polyshape| object can be easily plotted.
%
% *Note*: the duplicate points warning below can be safely ignored. If I
% were doing something other than a quick prototype, I would construct the
% data more carefully to avoid it.
%
% There you have it. You can plot a treble clef in MATLAB.
%
% What will you do with your new-found power?

p = polyshape(x,y);
plot(p)
axis equal
axis ij

%%
function placelabel(pt,str)
% Utility function for code demonstrating Bezier curve plotting
    x = pt(1);
    y = pt(2);
    h = line(x,y);
    h.Marker = '.';
    h = text(x,y,str);
    h.HorizontalAlignment = 'center';
    h.VerticalAlignment = 'bottom';
end

%%
% _Copyright 2020 The MathWorks, Inc._

% for i=1:length(clef_path)
%     
% end