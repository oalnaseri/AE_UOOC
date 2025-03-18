function [loss, gradientsEncoder, gradientsDecoder] = modelGradients(encoderNet, decoderNet, XBatch, YBatch, aU, bU, cU, h, SNR, Nt, Nr, FlipFlag, Eb_N0_dB)
    % Forward pass through the encoder
    encodedData = forward(encoderNet, XBatch);

    % Pass the encoded data through the UWOC channel
    channelOutput = real(UWOC_channel_ae(aU, bU, cU, extractdata(encodedData), h, SNR, Nt, Nr, FlipFlag, Eb_N0_dB));
    channelOutput = reshape(channelOutput, [size(encodedData, 1), size(encodedData, 2)]);
    channelOutput = dlarray(channelOutput, 'CB');

    % Forward pass through the decoder
    decodedData = forward(decoderNet, channelOutput);

    % Compute the loss (cross-entropy loss)
    loss = crossentropy(decodedData, YBatch);

    % Backward pass (compute gradients)
    gradientsEncoder = dlgradient(loss, encoderNet.Learnables);
    gradientsDecoder = dlgradient(loss, decoderNet.Learnables);
end