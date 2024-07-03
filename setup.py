from setuptools import setup

setup(
    name='PROMINENT_methylation',
    version='0.1.9',
    packages=['PROMINENT_methylation'],
    description='PROMINENT: AN INTERPRETABLE DEEP LEARNING METHOD TO PREDICT PHENOTYPES USING DNA METHYLATION',
    url='https://github.com/cloudmacchiato/dlmethylation',
    author='Laizhi Zhang',
    author_email='laz64@pitt.edu',
    license='MIT',
    entry_points={
        'console_scripts': [
            'PROMINENT-data_prepare = PROMINENT_methylation.dataprep:dataprep',
            'PROMINENT-train_test_cv = PROMINENT_methylation.train_test:train',
            'PROMINENT-scores = PROMINENT_methylation.scores:get_scores',
            'PROMINENT-model_interpret = PROMINENT_methylation.interpret:get_feature_importance',
            'PROMINENT-model_train_test_independent = PROMINENT_methylation.independent_test:independent_test'
        ]
    },
    install_requires=['numpy',
                      'shap',
                      'seaborn',
                      'matplotlib',
                      'torch',
                      'imblearn'
    ]
)
