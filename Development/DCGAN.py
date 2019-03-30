import keras
import numpy as np
import pandas as pd

def DCGAN_PSC(choice, filename):
    '''
    A convolutional based generative adversarial network that generates novel
    unnatural protein topologies. This network uses the phi, psi angles as
    well as the Ca distances as protein structure features
    '''
    CSTmax = 88.731
    shape = (150, 3)
    latent = 100
    batchs = 128
    epochs = 20000
    Disclr = 1e-7
    Advrlr = 1e-7
    data = pd.read_csv(filename, ';')
    phi = data[data.columns[2::3]].values
    psi = data[data.columns[3::3]].values
    cst = data[data.columns[4::3]].values
    phi /= 360
    psi /= 360
    cst /= CSTmax
    X = np.array([phi, psi, cst])
    X = np.swapaxes(X, 0, 2)
    X = np.swapaxes(X, 0, 1)
    D = keras.models.Sequential()
    D.add(keras.layers.Conv1D(64, kernel_size=3, input_shape=shape))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(64, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.MaxPooling1D())
    D.add(keras.layers.Conv1D(128, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(128, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.MaxPooling1D())
    D.add(keras.layers.Conv1D(256, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(256, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(256, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(256, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
#	D.add(keras.layers.MaxPooling1D())
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
#	D.add(keras.layers.MaxPooling1D())
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Conv1D(512, kernel_size=3))
    D.add(keras.layers.LeakyReLU(alpha=0.2))
    D.add(keras.layers.Flatten())
    D.add(keras.layers.Dense(1, activation='sigmoid'))
    G = keras.models.Sequential()
    G.add(keras.layers.Dense(79*3, activation='relu', input_dim=latent))
    G.add(keras.layers.Reshape((79, 3)))
    G.add(keras.layers.Conv1D(128, kernel_size=3))
    G.add(keras.layers.Activation('relu'))
    G.add(keras.layers.UpSampling1D())
    G.add(keras.layers.Conv1D(64, kernel_size=3))
    G.add(keras.layers.Activation('relu'))
    G.add(keras.layers.Conv1D(3, kernel_size=3))
    G.add(keras.layers.Activation('tanh'))
    DM = keras.models.Sequential()
    DM.add(D)
    DM.compile(optimizer=keras.optimizers.Adam(Disclr), loss='mean_squared_error', metrics=['accuracy'])
    AM = keras.models.Sequential()
    AM.add(G)
    AM.add(D)
    AM.compile(optimizer=keras.optimizers.Adam(Advrlr), loss='mean_squared_error', metrics=['accuracy'])
    if choice == 'train':
        for epoch in range(epochs):
            real = X[np.random.randint(0, X.shape[0], size=batchs)]
            noise = np.random.uniform(0.0, 1.0, size=[batchs, 100])
            fake = G.predict(noise)
            x = np.concatenate((real, fake))
            y = np.ones([2*batchs, 1])
            y[batchs:, :] = 0
            d_loss = DM.train_on_batch(x, y)
            y = np.ones([batchs, 1])
            a_loss = AM.train_on_batch(noise, y)
            D_loss = round(float(d_loss[0]), 3)
            D_accu = round(float(d_loss[1]), 3)
            A_loss = round(float(a_loss[0]), 3)
            if epoch % 500==0: print('Epoch {:5}\tDiscriminator loss: {:.7f}   |   Adversarial loss: {:.7f}'.format(epoch, D_loss, A_loss))
            G.save_weights('weights.h5')
        print('---------------------------------------------------------\nDone')
    elif choice == 'generate':
        G.load_weights('weights.h5')
        noise = np.random.normal(0.5, 0.5, (1, 100))
        gen = G.predict(noise)
        gen = gen.reshape([450])
        gen = np.ndarray.tolist(gen)
        phiout = gen[0::3]
        psiout = gen[1::3]
        cstout = gen[2::3]
        phiout = [x*360.0 for x in phiout]
        psiout = [x*360.0 for x in psiout]
        cstout = [x*CSTmax for x in cstout]
        savedata = open('prediction.txt', 'w')
        for p, s ,c in zip(phiout, psiout, cstout):
            savedata.write('{};{};{}\n'.format(p, s, c))
        savedata.close()

DCGAN_PSC('train', 'PSC_Helix_5.csv', 88.731)
DCGAN_PSC('generate', 'PSC_Helix_500.csv', 88.731)
