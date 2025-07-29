from setuptools import setup, find_packages

setup(
    name="cladeaid",
    version="0.1",
    description="LCA assignment and refining tools",
    author="Nicolas Locatelli",
    license="MIT",
    packages=find_packages(),  # finds the cladeaid/ package
    include_package_data=True,
    install_requires=[],
    entry_points={
        "console_scripts": [
            "cladeaid=cladeaid.__main__:main",
        ]
    },
)