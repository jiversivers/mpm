function switchButton(button, editfield)
if ~button.Value
    editfield.BackgroundColor = [1 1 1];
    editfield.Editable = 'on';
else
    editfield.BackgroundColor = [0.7 0.7 0.7];
    editfield.Editable = 'off';
end
end