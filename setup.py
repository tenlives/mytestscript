 # coding=utf-8
from setuptools import setup, find_packages
setup(
	name='error',
	version='0.1',
	packages = find_packages(),
	author='halazi',
	author_email='3033263880@qq.com',
	entry_points ={
		'console_scripts':[
			'error = demo.test:main',
		],
	},
	zip_safe=False
)
