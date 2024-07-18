% Signture experience Project - Air hockey Game
% Jeff Sheng

clc
clear

%% The rink layout
scrnsize = get(0,'ScreenSize');
x = 50;
y = 50;
w = scrnsize(3) - 2*x;
h = scrnsize(4) - 2.5*y;
fig = figure('Position', [x,y,w,h],'Color','cyan','Toolbar','none');
axl = axes('Position',[0.075, 0.075, 0.85, 0.85]);
axis equal;
axis(axl,[-1 1 -0.5 0.5]);
axis manual;
axis off
% Create rectangle with rouned corners and white center
rectangle('Position',[-1,-0.5,2,1],'cur',[0.05,0.1],'Linewidth',6,'Facecolor','w');
% Create red circle using rectangle function with curvarure set to [1,1]
rectangle('Position',[-0.125, -0.125, 0.25,0.25],'Curvature',[1,1],'LInewidth',3,'EdgeColor','red');
line([0,0],[-0.5,0.5],'Color','red','Linewidth',5,'LineStyle','-');
% Left goal
rectangle('Position',[-1.125, -0.15, 0.28,0.3],'Curvature',[1,1],'LInewidth',3,'EdgeColor','red');
% Right goal
rectangle('Position',[0.85, -0.15, 0.28,0.3],'Curvature',[1,1],'LInewidth',3,'EdgeColor','red');

%% Creat puck
% initalize the global variable
global quitgame blk_1_x blk_1_y blk_2_x blk_2_y ScoreBox_1 ScoreBox_2
kdl = @KeyDownListener;
set(fig, 'KeyPressFcn', @KeyDownListener)

% Create puck patch
r_puck = 0.075; % radius of puck
theta = linspace(0,2*pi);
xp_coor = r_puck * cos(theta);
yp_coor = r_puck * sin(theta);
puck = patch(xp_coor,yp_coor,'b');
pk_y = 0; %inital position of puck
% the puck would randomly be surfed by players
position = randi(2);
if position == 1
    pk_x = 0.5;
else
    pk_x = -0.5;
end
% puck initial velocity
pk_vx = 0;
pk_vy = 0;
dt = 0.02;
pk_mass = 1;
% limit of puck
pk_y_lim = 0.5 - r_puck;
pk_x_lim = 1 - r_puck;
%% Creat blockers
% blocker 1
blk_1_x = 0.7;% inital position of blocker 1
blk_1_y = 0;
r_block = 0.05;
xb_coor_1 = r_block * cos(theta);
yb_coor_1 = r_block * sin(theta);
block_1 = patch(xb_coor_1, yb_coor_1,'r');
% blocker 2
blk_2_x = -0.7;% inital postion of blocker 2
blk_2_y = 0;
xb_coor_2 = r_block * cos(theta);
yb_coor_2 = r_block * sin(theta);
block_2 = patch(xb_coor_2, yb_coor_2,'g');
% blocker volcity
blk_mass = 2; 
blk_vx = 0.1;
blk_vy = 0;
% limit of blocker
blk_y_lim = 0.5 - r_block;
blk_x_lim = 1 - r_block;
%% scoreborad
% Initialize score
score_1 = 0;
score_2 = 0;

% Create ScoreBox_1
ScoreBox_1 = uicontrol('Style', 'edit', 'String',score_1, ...
    'Position', [500, 870, 100, 20]);
% Create label on Player 1 score
uicontrol('Style', 'text', 'String', 'Player 1 ', ...
    'Position', [500, 890, 100, 20])
% Create ScoreBox_2
ScoreBox_2 = uicontrol('Style', 'edit', 'String',score_2, ...
    'Position', [1200, 870, 100, 20]);
% Create label on Player 2 score
uicontrol('Style', 'text', 'String', 'Player 2 ', ...
    'Position', [1200, 890, 100, 20])
%% the core movement
% Initialize variable quitgame
quitgame = 0;
while quitgame ~= 1
    % set the limit of the puck and blocker
    [pk_vx, pk_vy, pk_x, pk_y,score_1, score_2] = puck_limit(pk_x, pk_y,...
        pk_vx, pk_vy, pk_x_lim, pk_y_lim,score_1, score_2)
    % call the block limit function
    block_limit(blk_x_lim, blk_y_lim, r_block)
    % call the collision funcion on blocker 1
    dis_blk_1 = sqrt((blk_1_x - pk_x)^2 + (blk_1_y - pk_y)^2);
    if dis_blk_1 < r_block + r_puck 
        [pk_x,pk_y,pk_vy,pk_vx] = Collision(blk_1_x,blk_1_y,blk_vx,blk_vy,...
        r_block,blk_mass,pk_x,pk_y,pk_vx,pk_vy,r_puck,pk_mass)
    end
    % call the collision function on blocker 2
    dis_blk_2 = sqrt((blk_2_x - pk_x)^2 + (blk_2_y - pk_y)^2);
    if dis_blk_2 < r_block + r_puck 
        [pk_x,pk_y,pk_vy,pk_vx] = Collision(blk_2_x,blk_2_y,blk_vx,blk_vy,...
        r_block,blk_mass,pk_x,pk_y,pk_vx,pk_vy,r_puck,pk_mass)
    end
    % puck movement
    pk_x = pk_x + pk_vx * dt;
    pk_y = pk_y + pk_vy * dt;
    set(puck,'XData', xp_coor + pk_x, 'YData', yp_coor + pk_y);
    % use the keypress function to control the blockers
    set(block_1, 'XData', xb_coor_1 + blk_1_x, 'YData', yb_coor_1 + blk_1_y);
    set(block_2, 'XData', xb_coor_2 + blk_2_x, 'YData', yb_coor_2 + blk_2_y);
    pause(0.01);
    % the score up to 7 
    if score_1 == 7 
        quitgame = 1;
    end
    if score_2 == 7
        quitgame = 1;
    end
    
end   
close all
clear all


%% Keypress Function
function KeyDownListener(src, event)
% KeyDownListener is triggered by pressing a key on the keyboard
% Inputs:  src   - the source of the keypress event
%          event - the actual keypress event
% These inputs are required exactly as shown
% No output variables - must be global variables

% Give function access to global variables for quit 'flag' & position
global quitgame blk_1_x blk_1_y blk_2_x blk_2_y 


% Read name of key from key-press event that triggered this function
keyID = event.Key;

% Try switch-case structure to set global 'output' variable(s)
switch keyID 
    case 'q'
        quitgame = true;
    case 'uparrow'
        % Move block 1 up by 1 unit
        blk_1_y = blk_1_y + 0.025;
    case 'downarrow'
        % Move block 1 down by 1 unit
        blk_1_y = blk_1_y - 0.025;
    case 'rightarrow'
        % Move block 1 right by one unit
        blk_1_x = blk_1_x + 0.025;
    case 'leftarrow'
        % Move block 1 left by one unit
        blk_1_x = blk_1_x - 0.025;
    case 'w' 
        % Move blokc 2 up by 1 unit
        blk_2_y = blk_2_y + 0.025;
    case 's'
        % Move block 2 down by 1 unit
        blk_2_y = blk_2_y - 0.025;
    case 'd' 
        % Move block 2 right by 1 unit
        blk_2_x = blk_2_x + 0.025;
    case 'a'
        % Move block 2 left by 1 unit
        blk_2_x = blk_2_x - 0.025;
        
    otherwise
        keyID
end

end      

%% puck limit function
function [pk_vx, pk_vy, pk_x, pk_y,score_1, score_2] = puck_limit...
    (pk_x, pk_y,pk_vx, pk_vy, pk_x_lim, pk_y_lim,score_1, score_2) 
global ScoreBox_1 ScoreBox_2
% set the top limit of puck
if pk_y > pk_y_lim 
    pk_vy = -pk_vy;
    pk_y = pk_y_lim;
end
% set the button limit of puck
if pk_y < -pk_y_lim
    pk_vy = -pk_vy;
    pk_y = -pk_y_lim;
end
% puck limit on the right side
if pk_x > pk_x_lim 
    % set the goal for player 1
    if pk_y > -0.15 && pk_y < 0.15 
        score_1 = score_1 + 1;
        ScoreBox_1.String = score_1;
        % reset the puck 
        pk_x = 0.5;
        pk_y = 0;
        pk_vx = 0;
        pk_vy = 0;
    else 
        % puck limit
        pk_vx = -pk_vx;
        pk_x = pk_x_lim;
    end
end
% set the puck limit on the left side
if pk_x < -pk_x_lim
    % set the goal for player 2
    if pk_y > -0.15 && pk_y < 0.15
        score_2 = score_2 + 1;
        ScoreBox_2.String = score_2;
        % reset the puck
        pk_x = -0.5;
        pk_y = 0;
        pk_vx = 0;
        pk_vy = 0;
    else 
        % puck limit
        pk_vx = -pk_vx;
        pk_x = -pk_x_lim;
    end
end

end
%% block_limit function
function block_limit(blk_x_lim, blk_y_lim, r_block)
    global blk_1_x blk_1_y blk_2_x blk_2_y
    % the left limit of blocker 1
    if blk_1_x < r_block
       blk_1_x = r_block;
    end
    % the right limit of blocker 1
    if blk_1_x > blk_x_lim
       blk_1_x = blk_x_lim;
    end
    % the top limit of blocker 1
    if blk_1_y > blk_y_lim
       blk_1_y = blk_y_lim;
    end
    % the button limit of blocker 1
    if blk_1_y < -blk_y_lim
       blk_1_y = -blk_y_lim;
    end
    % the right limit of blocker 2
    if blk_2_x > -r_block
       blk_2_x = -r_block;
    end
    % the left limit of blocker 2
    if blk_2_x < -blk_x_lim
       blk_2_x = -blk_x_lim;
    end
    % the top limit of blocker 2
    if blk_2_y > blk_y_lim
       blk_2_y = blk_y_lim;
    end
    % the button limit of blocker 2
    if blk_2_y < -blk_y_lim
       blk_2_y = -blk_y_lim;
    end
end   

%% Collision function
% Function to calculate the velocity components of the puck after a
% collision with a blocker in an air hockey game.
% Inputs:
%         xb  --   x position of the blocker
%         yb  --   y position of the blocker
%         Vxb --   x velocity of the blocker
%         Vyb --   y velocity of the blocker
%         rb  --   radius of the blocker
%         mb  --   mass of the blocker
%         xp  --   x position of the puck
%         yp  --   y position of the puck
%         Vxp --   x velocity of the puck
%         Vyp --   y velocity of the puck
%         rp  --   radius of the puck
%         mp  --   mass of the puck

function [xp,yp,Vyp,Vxp] = Collision(xb,yb,Vxb,Vyb,rb,mb,xp,yp,Vxp,Vyp,rp,mp)
    % Calculate the angle between the x-axis and the line between the
    % center of the puck and center of the blocker. The axis between the
    % two centers is the normal axis (normal to the collision) and the
    % tangential axis is perpendicular to the normal axis at the contact
    % point between the blocker and the puck.  The positive direction of
    % the normal axis is from the puck to the blocker.
%     xpin=xp;                  %These values are here to help debug
%     ypin=yp;                  %because the values get calculated as 
%     Vxpin = Vxp;              %outputs later in the code
%     Vypin = Vyp;
    theta = -atan2d((yb-yp),(xb-xp));
    % Calculate the angle between the x-axis and the velocity vector  
    %    of the blocker.        
    % alpha = atan2d(Vyb,Vxb);
    alpha = 180 - theta;
    % Calculate the angle between the x-axis and the velocity vector  
    %    of the puck.        
    beta = atan2d(Vyp,Vxp);
    % Move the puck along the normal axis in the negative direction to
    % eliminate any overlap between the puck and the blocker
    xp = xb - (rp+rb) * cosd(theta);        
    yp = yb + (rp+rb) * sind(theta);
    % Calculate the total velocity of the blocker and the puck
    Vb1 = sqrt(Vyb^2 + Vxb^2);
    Vp1 = sqrt(Vyp^2 + Vxp^2);
    % Calculate the components of the blocker velocity in the normal and
    % tangential directions ( n and s dirctions )
    Vb1n = Vb1 * cosd(theta + alpha);
    Vb1s = Vb1 * sind(theta + alpha);
    % Calculate the components of the puck velocity in the normal and
    % tangential directions ( n and s dirctions ) before the collision
    Vp1n = Vp1 * cosd(theta + beta);
    Vp1s = Vp1 * sind(theta + beta);
    % Since there is no friction in our case, there can be no change in
    % velocity in the tangential direction, thus the tangential velocity of
    % the blocker and puck are:
    Vb2s = Vb1s;
    Vp2s = Vp1s;
    % Now we must calculate the velocity of the puck and the blocker in the
    % normal direction after the collision.  To do this we must balance the
    % the total momentum before the collision with the total momentum after
    % the collision.  Since we are saying that the collision is elastic,
    % the kinetic energy after the collision must also equal the kinetic
    % energy before the collision
    % Calculate the total momentum before the collision
    P1n = mp * Vp1n + mb * Vb1n; 
    % Calculate the total momentum after the collision
    KE1 = 0.5*mp*Vp1^2 + 0.5*mb*Vb1^2;
    %
    % The kinetic energy after the collision is
    %            KE2 = 0.5*mp*Vp2^2 + 0.5*mb*Vb2^2;
    %            KE1 = 0.5*mp*Vp2^2 + 0.5*mb*Vb2^2;
    %          2*KE1 = mp*Vp2^2 + mb*Vb2^2;
    % Substituting from above for Vb2 and rearranging terms we can put this
    % equation in the following form a*Vp2n^2 + b * Vp2n + c = 0 where the
    % coefficients a, b, and c are shown below: 
    a = (mp^2 + mp*mb)/mb;
    b = -(2*P1n*mp)/mb;
    c = (P1n^2/mb) + mp*Vp1s^2 + mb*Vb1s^2 - 2*KE1;
    % We can then solve for Vp2n using the quadratic formula. (Note that
    % since the two bodies contact along the normal axis, the Vp2n must be
    % more negative than Vp1n, thus we select the negative value for the
    % square root.
    Vp2n = (-b-sqrt(b^2 - 4*a*c))/(2*a);
    % The value for Vb2n can be found from substituting back into the
    % equation for the momentum balance, though this is actually
    % unnecessary, since we do not need to calculate the velocity of the
    % blocker.  
    Vb2n = (P1n - mp*Vp2n)/mb;
    % The magnitude of the velocity of the puck can be found by taking the
    % square root of the sum of the squares of the components in the normal
    % and tangential directions.
    Vp2 = sqrt(Vp2n^2 + Vp2s^2);
    % The new angle for the puck can be found by calculating the angle with
    % respect to the normal and tangential axes and then subtracting the 
    % angle between the x and y axes and the normal and tangential axes
    beta2 = atan2d(Vp2s,Vp2n) - theta;
    % The magnitude of the velocity of the blocker can be found by taking
    % the square root of the sum of the squares of the components in the
    % normal and tangential directions.
    Vb2 = sqrt(Vb2n^2 + Vb2s^2);
    % The new angle for the blocker can be found by calculating the angle 
    % with respect to the normal and tangential axes and then subtracting  
    % the angle between the x and y axes and the normal and tangential axes
    alpha2 = atan2d(Vb2s,Vb2n) - theta;
    % The new velocities in the x and y directions for the puck are just
    % the new puck velocity times the cos and sin of the angle beta2.
    Vyp = Vp2 * sind(beta2);
    Vxp = Vp2 * cosd(beta2);
end

    
    
    
