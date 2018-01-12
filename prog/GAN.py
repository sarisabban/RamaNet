import numpy , pandas , keras
#https://arxiv.org/abs/1511.06434 https://github.com/kyloon/dcgan

#Import data
data = pandas.read_csv('dataCA.csv' , ';')

#Convert data into a 3D numpy tensor, the shape will be (XYZ , number of examples , number of amino acids)
x = pandas.DataFrame.as_matrix(data[data.columns[2:]][data[data.columns[2:]].columns[::3]])
y = pandas.DataFrame.as_matrix(data[data.columns[3:]][data[data.columns[3:]].columns[::3]])
z = pandas.DataFrame.as_matrix(data[data.columns[4:]][data[data.columns[4:]].columns[::3]])
xyz = numpy.array([x , y , z])
print(xyz.shape)						#(3, 20, 150)

def Generative():
	model = keras.models.Sequential()
#	model.add(keras.layers.Dense())
#	model.add(keras.layers.Dense())
#	model.add(keras.layers.BatchNormalization())
#	model.add(keras.layers.Reshape())
#	model.add(keras.layers.UpSampling1D())
#	model.add(keras.layers.Conv1D())
#	model.add(keras.layers.UpSampling1D())
#	model.add(keras.layers.Conv1D())
#	model.compile(keras.optimizers.Adam() , loss = 'binary_crossentropy' , metrics = ['accuracy'])
	return(model)

#Generative().summary()

def Discriminative():
	model = keras.models.Sequential()
	model.add(keras.layers.Conv1D(1 , 5 , input_shape = (150 , 3) , activation = 'relu'))
	model.add(keras.layers.MaxPooling1D((2)))
	model.add(keras.layers.Conv1D(1 , 5 , activation = 'relu'))
	model.add(keras.layers.MaxPooling1D((2)))
	model.add(keras.layers.core.Flatten())
	model.add(keras.layers.Dense(34 , activation = 'relu'))
	model.add(keras.layers.Dense(2 , activation = 'softmax'))
	model.compile(keras.optimizers.Adam() , loss = 'binary_crossentropy' , metrics = ['accuracy'])
	return(model)

Discriminative().summary()

def GAN(G , D):
	model = keras.models.Sequential()
	model.add(G)
	D.trainable = False
	model.add(D)
#	model.compile(keras.optimizers.Adam() , loss = 'binary_crossentropy' , metrics = ['accuracy'])
	return(model)

#GAN(Generative() , Discriminative())

'''
GAN To Design Molecules
I am looking for a freelancer to help me (code and advise) develop a Generative Adversarial Network (GAN) that generates spatial atom positions (XYZ) to draw molecules.
I have generated the nessesary dataset and I have part of the code roughly written on Keras.
This is part of a research project, so you can choose to either have payment for your services, or have your name as a co-author in the resultant published paper.
I have access to a High Preformance Cluster computer to train the model, so you will not need to setup or use any cloud computing service.
'''



'''
def train(batch_size, num_epoch):
    (X_train, y_train), (X_test, y_test) = mnist.load_data()
    X_train = (X_train.astype(np.float32) - 127.5)/127.5
    X_train = X_train.reshape((X_train.shape[0], 1) + X_train.shape[1:])
    discriminator = discriminator_model()
    generator = generator_model()
    discriminator_on_generator = \
        generator_containing_discriminator(generator, discriminator)
    d_optim = SGD(lr=0.0005, momentum=0.9, nesterov=True)
    g_optim = SGD(lr=0.0005, momentum=0.9, nesterov=True)
    generator.compile(loss='binary_crossentropy', optimizer="SGD")
    discriminator_on_generator.compile(
        loss='binary_crossentropy', optimizer=g_optim)
    discriminator.trainable = True
    discriminator.compile(loss='binary_crossentropy', optimizer=d_optim)
    noise = np.zeros((batch_size, 100))
    for epoch in range(num_epoch):
        print("Epoch " + str(epoch+1) + "/" + str(num_epoch) +" :")
        print("Number of batches:", int(X_train.shape[0]/batch_size))
        for index in range(int(X_train.shape[0]/batch_size)):
            for i in range(batch_size):
                noise[i, :] = np.random.uniform(-1, 1, 100)
            image_batch = X_train[index*batch_size:(index+1)*batch_size]
            generated_images = generator.predict(noise, verbose=0)
            if index % 20 == 0:
                image = combine_images(generated_images)
                image = image*127.5+127.5
                Image.fromarray(image.astype(np.uint8)).save(
                    "images/"+str(epoch+1)+"_"+str(index+1)+".png")
            X = np.concatenate((image_batch, generated_images))
            y = [1] * batch_size + [0] * batch_size
            d_loss = discriminator.train_on_batch(X, y)
            print("Batch %d d_loss : %f" % (index+1, d_loss))
            for i in range(batch_size):
                noise[i, :] = np.random.uniform(-1, 1, 100)
            discriminator.trainable = False
            g_loss = discriminator_on_generator.train_on_batch(
                noise, [1] * batch_size)
            discriminator.trainable = True
            print("Batch %d g_loss : %f" % (index+1, g_loss))
            if index % 10 == 9:
                generator.save_weights('generator_weights', True)
                discriminator.save_weights('discriminator_weights', True)

def generate(batch_size, pretty=False):
    generator = generator_model()
    generator.compile(loss='binary_crossentropy', optimizer="SGD")
    generator.load_weights('generator_weights')
    if pretty:
        discriminator = discriminator_model()
        discriminator.compile(loss='binary_crossentropy', optimizer="SGD")
        discriminator.load_weights('discriminator_weights')
        noise = np.zeros((batch_size*20, 100))
        for i in range(batch_size*20):
            noise[i, :] = np.random.uniform(-1, 1, 100)
        generated_images = generator.predict(noise, verbose=1)
        d_pret = discriminator.predict(generated_images, verbose=1)
        index = np.arange(0, batch_size*20)
        index.resize((batch_size*20, 1))
        pre_with_index = list(np.append(d_pret, index, axis=1))
        pre_with_index.sort(key=lambda x: x[0], reverse=True)
        pretty_images = np.zeros((batch_size, 1) +
                               (generated_images.shape[2:]), dtype=np.float32)
        for i in range(int(batch_size)):
            idx = int(pre_with_index[i][1])
            pretty_images[i, 0, :, :] = generated_images[idx, 0, :, :]
        image = combine_images(pretty_images)
    else:
        noise = np.zeros((batch_size, 100))
        for i in range(batch_size):
            noise[i, :] = np.random.uniform(-1, 1, 100)
        generated_images = generator.predict(noise, verbose=1)
        image = combine_images(generated_images)
    image = image*127.5+127.5
    Image.fromarray(image.astype(np.uint8)).save(
"images/generated_image.png")
'''
