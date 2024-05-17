from setuptools import setup

setup(
    name='prominent',
    version='0.1',
    packages=['prominent'],
    description='PROMINENT: AN INTERPRETABLE DEEP LEARNING METHOD TO PREDICT PHENOTYPES USING DNA METHYLATION',
    url='https://github.com/cloudmacchiato/dlmethylation',
    author='Laizhi Zhang',
    author_email='laz64@pitt.edu',
    license='MIT',
    entry_points={
        'console_scripts': [
            'prominent-data_prepare = prominent.dataprep:dataprep',
            'prominent-train_test = prominent.train_test:train',
            'prominent-scores = prominent.scores:get_scores',
            'prominent-model_interpret = prominent.interpret:get_feature_importance'
        ]
    },
    install_requires=['numpy',
                      'pandas',
                      'scikit-learn',
                      'shap',
                      'pickle',
                      'seaborn',
                      'matplotlib',
                      'torch',
                      'imblearn'
    ]
)