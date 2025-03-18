clc;
clear;
close all;

% Define the dimensions
inputSize = 1; % Size of the input (1x10000 array)
hiddenSize1 = 64; % Number of units in the first hidden layer
hiddenSize2 = 1; % Number of units in the second hidden layer (bottleneck)
outputSize = 2; % Number of classes (binary classification)

% Encoder network
encoderLayers = [
    featureInputLayer(inputSize, 'Name', 'input') % Input layer (no sequence dimension)
    fullyConnectedLayer(hiddenSize1, 'Name', 'enc_fc1', ...
        'Weights', randn(hiddenSize1, inputSize) * 0.01, ... % Initialize weights
        'Bias', zeros(hiddenSize1, 1)) % Initialize biases
    eluLayer('Name', 'enc_elu1') % ELU activation
    fullyConnectedLayer(hiddenSize2, 'Name', 'enc_fc2', ...
        'Weights', randn(hiddenSize2, hiddenSize1) * 0.01, ... % Initialize weights
        'Bias', zeros(hiddenSize2, 1)) % Initialize biases
];
encoderNet = dlnetwork(encoderLayers);

% Decoder network
decoderLayers = [
    featureInputLayer(hiddenSize2, 'Name', 'dec_input') % Input layer for decoder
    fullyConnectedLayer(hiddenSize1, 'Name', 'dec_fc1', ...
        'Weights', randn(hiddenSize1, hiddenSize2) * 0.01, ... % Initialize weights
        'Bias', zeros(hiddenSize1, 1)) % Initialize biases
    eluLayer('Name', 'dec_elu1') % ELU activation
    fullyConnectedLayer(outputSize, 'Name', 'dec_fc2', ...
        'Weights', randn(outputSize, hiddenSize1) * 0.01, ... % Initialize weights
        'Bias', zeros(outputSize, 1)) % Initialize biases
    softmaxLayer('Name', 'dec_softmax') % Softmax layer for classification
];
decoderNet = dlnetwork(decoderLayers);

% Generate some random input data for demonstration
numSamples = 10000;
X = 17.9310195815304 * (2 * randi([0, 1], 1, numSamples) - 1); % Randomly generate -17.93 or 17.93
X = reshape(X, [inputSize, numSamples]); % Reshape to [numFeatures x numSamples]
X = dlarray(X, 'CB'); % Convert to dlarray (1x10000)

% Convert input data to categorical labels
Y = categorical(X == 17.9310195815304, [false, true], {'Class0', 'Class1'}); % Map to categorical labels
Y = onehotencode(Y, 1); % Convert to one-hot encoded labels

% Define UWOC channel parameters
aU = 0.1; % Absorption coefficient
bU = 0.2; % Scattering coefficient
cU = 0.3; % Attenuation coefficient
h = 1; % Channel impulse response (simplified for this example)
SNR = 10; % Signal-to-noise ratio
Nt = 1; % Number of transmit antennas
Nr = 1; % Number of receive antennas
FlipFlag = 0; % Flag for flipping the signal
Eb_N0_dB = 10; % Energy per bit to noise power spectral density ratio

% Training options
numEpochs = 100;
miniBatchSize = 32;
learningRate = 0.001;

% Training loop
for epoch = 1:numEpochs
    fprintf('Epoch %d/%d\n', epoch, numEpochs);
    
    % Shuffle data
    idx = randperm(numSamples);
    X = X(:, idx);
    Y = Y(:, idx);
    
    % Mini-batch training
    for i = 1:miniBatchSize:numSamples
        % Get mini-batch
        XBatch = X(:, i:min(i+miniBatchSize-1, numSamples));
        YBatch = Y(:, i:min(i+miniBatchSize-1, numSamples));

        % Compute gradients and loss using dlfeval
        [loss, gradientsEncoder, gradientsDecoder] = dlfeval(@modelGradients, encoderNet, decoderNet, XBatch, YBatch, aU, bU, cU, h, SNR, Nt, Nr, FlipFlag, Eb_N0_dB);
        
        % Update network parameters
        encoderNet.Learnables = dlupdate(@(params, grad) params - learningRate * grad, encoderNet.Learnables, gradientsEncoder);
        decoderNet.Learnables = dlupdate(@(params, grad) params - learningRate * grad, decoderNet.Learnables, gradientsDecoder);
        
        % Display loss
        fprintf('Loss: %.4f\n', extractdata(loss));
    end
end

% Test the autoencoder with some data
testData = X(:, 1:5); % Use the first 5 samples for testing
testLabels = Y(:, 1:5); % Corresponding labels

% Get encoded data (bottleneck representation)
encodedData = forward(encoderNet, testData);

% Pass the encoded data through the UWOC channel
channelOutput = real(UWOC_channel_ae(aU, bU, cU, extractdata(encodedData), h, SNR, Nt, Nr, FlipFlag, Eb_N0_dB));
channelOutput = reshape(channelOutput, [hiddenSize2, size(encodedData, 2)]);
channelOutput = dlarray(channelOutput, 'CB');

% Get decoded data (reconstructed output)
decodedData = forward(decoderNet, channelOutput);

% Display the results
disp('Original Data:');
disp(extractdata(testData));
disp('Original Labels:');
disp(testLabels');
disp('Decoded Data:');
disp(extractdata(decodedData));

% Function to compute gradients
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