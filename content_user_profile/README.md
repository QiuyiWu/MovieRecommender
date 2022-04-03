## Content-based Movie Recommendation System using User Profiles
----

### Problem Formulation and Data Pipeline

Here we have three entities, user $U$, movie $M$, and rating $R$. Any user can rate any arbitrary number of movies; therefore, the amount of ratings $|\{R\}|$ we have is $|\{U \times M\}|$ which is around 25M.

We constructed a movie profile based on the given 1,128 tag genome scores, $d=1128$, and some preprocessing. Namely, we write a movie profile $\vec{M} \in \mathbb{R}^d$ with TF-IDF transformation as below.
$$ \vec{M_i} = tfidf(M_i, \{M\}) $$

And the user profile is constructed as follow.
$$ \vec{U}(U, \{M_j | U\}) = \frac{\sum_{j}^{|\{M_j | U\}|} {R(U, M_j) \cdot \vec{M_j}}} {\sum_{j}^{|\{M_j | U\}|} {R(U, M_j)}}$$
Despite the ugly notation above, in essence, every user profile is a normalized weighted sum of all the movies that a user had rated in the pasted. The weights are chosen as the ratings of this user gave to all the movies.

Lastly, the model is trying to predict the new rating given a user profile and a movie profile,
$$ <\vec{U}, \vec{M}> \Rightarrow R(U, M) $$
where $ <\vec{U}, \vec{M}> $ is a feature vector that is $\in \mathbb{R}^{2d}$.

---

### Create the virtual enviornment using Conda

First, you will need to install miniconda to use it for managing your environments. [Here is the official Anaconda documentation](https://docs.conda.io/en/latest/miniconda.html).

Now go to your favorite CLI, move to this folder directory, and run the below command to create the conda environment locally under the folder `.env`. Furthermore, the `environment.yml` is the YAML file containing all python package requirements you will need to run this example.
```bash
conda env create -p .env -f environment.yml
```

To activate the environment, do below. Now, you are done and ready to run the notebook
```bash
conda activate ./.env
```

----

### Where to get the data and pre-trained model

Next, we will need to get all the data to construct our training and testing set. The data was preprocessed based on MovieLens 25M, which you [can find it here](https://grouplens.org/datasets/movielens/). However, you do not need to download the original raw data. Here is [the google drive link that hosts the data set](https://drive.google.com/drive/folders/1F7kHWFI8KMOe_c_vVn462zBYH0jXqd6r?usp=sharing). We will be using our Linear Model with Stochastic Gradient Descent (SGD) in training.

Be sure to download the files from google drive to the `./data` folder in this directory. The notebook will be reading them from this specific directory.

----

### Data doc

We have, in total, three `.feather` binary files. The final training and testing dataset must be constructed based on these three files for three reasons.
- We are using SGD mini-bath methodology to train the linear model. Thus, we do it iteratively instead of feeding the entire dataset, which is about [25M x 2k] size.
- This entire dataset is too big to process at once for a typical personal computer.
- We can still achieve good results without the need to see all available data.

Here are the files and their documentation.
You will need to use pandas and pyarrow to read them like `pd.read_feather(file_name)`
- `feature_hist.feather`: it contains all the rating data from MovieLens 25M.
    - The `is_train` label is used to enforce reproducibility and prevent test data from creeping into the training set.
    - In addition, here, we only select movie ratings with tag genome features.
    - Columns: [userId: Int, movieId: Int, rating: Int, is_train: Boolean]
- `feature_movies.feather`: here, we have genome features for each movie. Or you can think of them as movie profiles.
    - Tag genome is a dense matrix that is 100% filled for each movie and consists of 1,128 different tags.
    - We preprocessed the genome score with TF-IDF and scaled it up by 100x to prevent precision error before constructing this dataset.
    - Columns [movieId: Int, genome: Array[Float]]
- `feature_profile.feather`: this is the user profile dataset
    - User profile is constructed based on movies profiles for each specific user.
    To create a train/test set, we have to do this separately to prevent the model from seeing test data during training. Because the user profile is an aggregation of all the movies a user has seen so far, we should not include some movies into the training set because they have not seen them yet.
    - User profile feature array should be the same size as any movie profile feature array. In this case, it is a vector with size $\mathbb{R}^{1128}$
    - Columns [userId: Int, train: Array[Float], test: Array[Float]]
- `model_sgd_regressor.joblib`: this is the `joblib` binary of the model from `sklearn`
- `model_sgd_regressor_scores.feather`: and here contains all the training and RMSE history
