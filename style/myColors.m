% Return Hex Code for colorblind safe combinations of colors

% 'Vibrant' Color Package from
% https://personal.sron.nl/~pault/#sec:qualitative

function [HEX] = myColors(color)
    
switch color
    case 'blue'
        HEX = '#0077BB';
    case 'cyan'
        HEX = '#33BBEE';
    case 'green'
        HEX = '#009988';
    case 'orange'
        HEX = '#EE7733';
    case 'red'
        HEX = '#CC3311';
    case 'magenta'
        HEX = '#EE3377';
    case 'grey'
        HEX = '#BBBBBB';
    case 'black'
        HEX = '#000000';
end
        
end