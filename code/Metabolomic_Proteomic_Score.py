import numpy as np
import pandas as pd
import keras
from keras.models import Sequential
from sklearn.model_selection import cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler  

from keras.layers import Dense, Dropout, Flatten, Conv1D, MaxPooling1D, AveragePooling1D, GRU, GlobalMaxPooling1D, GlobalAveragePooling1D
import matplotlib.pyplot as plt
# Import seaborn
import seaborn as sns
import tensorflow as tf

import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"

data = pd.read_csv('data/data_meta.csv')

data_cleaned = data.dropna(subset=['Heel_BMD_right_i0_T'])

X = data_cleaned.iloc[:, 16:]
y = data_cleaned['Heel_BMD_right_i0_T']

X_train, X_test, Y_train, Y_test = train_test_split(
    X, y, test_size=0.7, random_state=42
)

from sklearn.preprocessing import StandardScaler

scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv1D, MaxPooling1D, Flatten, Dense, Dropout, BatchNormalization
from tensorflow.keras import regularizers

model = Sequential()

# Block 1
model.add(Conv1D(32, 2, activation='selu',
                 input_shape=(X_train_scaled.shape[1], 1)))
model.add(BatchNormalization())
model.add(MaxPooling1D(2))

# Block 2
model.add(Conv1D(32, 2, activation='selu'))
model.add(BatchNormalization())
model.add(MaxPooling1D(2))

# Block 3
model.add(Conv1D(32, 2, activation='selu'))
model.add(MaxPooling1D(2))

model.add(Flatten())

# Dense layer
model.add(Dense(32, activation='selu',
                kernel_regularizer=regularizers.l2(1e-4)))
model.add(Dropout(0.3))

model.add(Dense(16, activation='selu',
                kernel_regularizer=regularizers.l2(1e-4)))

model.add(Dense(1, activation='linear'))

model.summary()

from tensorflow.keras.optimizers import Adam

optimizer = Adam(learning_rate=1e-4)

model.compile(
    optimizer=optimizer,
    loss='mse',
    metrics=['mae']
)


callback = keras.callbacks.EarlyStopping(monitor='loss', patience=10)
history = model.fit(
    X_train_scaled, Y_train,
    validation_split=0.2,
    epochs=100,
    batch_size=64,
    callbacks=[callback],
    shuffle=True   
)



plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('Model Loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend(['Train', 'Validation'])
plt.show()


Y_train_pred = model.predict(X_train_scaled)
Y_test_pred = model.predict(X_test_scaled)

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, r2_score

y_train_true = np.array(Y_train).flatten()
y_train_pred = np.array(Y_train_pred).flatten()

y_test_true = np.array(Y_test).flatten()
y_test_pred = np.array(Y_test_pred).flatten()

def calc_metrics(y_true, y_pred):
    r, _ = pearsonr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    return r, r2, rmse

r_train, r2_train, rmse_train = calc_metrics(y_train_true, y_train_pred)
r_test, r2_test, rmse_test = calc_metrics(y_test_true, y_test_pred)


plt.figure(figsize=(7,7))

# training
sns.scatterplot(
    x=y_train_true,
    y=y_train_pred,
    s=8,
    color='blue',
    alpha=0.5,
    label='Train'
)

# testing
sns.scatterplot(
    x=y_test_true,
    y=y_test_pred,
    s=8,
    color='black',
    alpha=0.6,
    label='Test'
)


min_val = min(y_train_true.min(), y_test_true.min(),
              y_train_pred.min(), y_test_pred.min())
max_val = max(y_train_true.max(), y_test_true.max(),
              y_train_pred.max(), y_test_pred.max())

plt.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=1)


textstr = '\n'.join((
    f'Train: r = {r_train:.2f}, R² = {r2_train:.2f}, RMSE = {rmse_train:.2f}',
    f'Test:   r = {r_test:.2f}, R² = {r2_test:.2f}, RMSE = {rmse_test:.2f}'
))

plt.text(
    0.05, 0.95,
    textstr,
    transform=plt.gca().transAxes,
    fontsize=11,
    verticalalignment='top',
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8)
)

plt.xlabel('True value', fontsize=12)
plt.ylabel('Predicted value', fontsize=12)
plt.title('Model performance (Train vs Test)', fontsize=13)
plt.legend()
plt.grid(alpha=0.2)

plt.tight_layout()
plt.show()


eid = data.iloc[:, 0]

y_train_true = np.array(Y_train).flatten()
y_train_pred = np.array(Y_train_pred).flatten()

y_test_true = np.array(Y_test).flatten()
y_test_pred = np.array(Y_test_pred).flatten()

train_eid = eid.loc[X_train.index]
test_eid = eid.loc[X_test.index]

df_train = pd.DataFrame({
    'eid': train_eid.values,
    'True': y_train_true,
    'Predicted': y_train_pred,
    'Dataset': 'Train'
})

df_test = pd.DataFrame({
    'eid': test_eid.values,
    'True': y_test_true,
    'Predicted': y_test_pred,
    'Dataset': 'Test'
})

df_all = pd.concat([df_train, df_test], axis=0)
df_all['Residual'] = df_all['True'] - df_all['Predicted']
