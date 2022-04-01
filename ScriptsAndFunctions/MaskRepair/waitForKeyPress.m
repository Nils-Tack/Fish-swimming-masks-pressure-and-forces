function keyChar = waitForKeyPress()
% Wait for a key press. Return the key pressed.

set(gcf, 'Windowkeypressfcn', @(src, event) set(src, 'UserData', true) );
set(gcf, 'UserData', false);     %reset key-has-been-pressed
while 1
    drawnow();    %give the graphics queue a chance to process the callback
    if get(gcf, 'UserData')   %was key pressed?
        hfig = gcf;
        keyChar = hfig.CurrentCharacter;
%         sprintf('Pressed: %s', keyChar)
        break;
    end
   pause(0.1) 
end