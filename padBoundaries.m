function B = padBoundaries(B)
% Pads boundaries to avoid NaN in iterpolation
B = cat(1, B(end,:,:), B, B(1,:,:));
B = cat(2, B(:,end,:), B, B(:,1,:));
B = cat(3, B(:,:,end), B, B(:,:,1));
end