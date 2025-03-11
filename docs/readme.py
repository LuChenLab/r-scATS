****** 快速上手
pip install sphinx

sphinx-quickstart# 快速开始
# 这里选择yes，输入y并回车。分离的好处是工程结构更加清晰。
# 项目名称
# 作者名字
# 版本
# 项目语种：en英文 ，zh_CN表示中文
### 这样就自动在当前目录生成基本文件：_build是存放输出的东西，_static是css相关，_templates是模版 渲染文件

make html # 输出网页版，在build下面生成相关文件 ， _build/html/index.html是我们的网页
# 每次影响网页结果的修改都需要重新运行这个命令，要不然网页也不会更新


****** markdown格式
# Sphinx默认支持reStructuredText，但是大部分人都是用Markdown来写文章，所以这里就需要安装支持Markdown的插件MyST-Parser
pip install -U myst-parser

## 详细参见：https://blog.csdn.net/whahu1989/article/details/140307613

****** 其他插件
# 1. 模版
pip install sphinx_rtd_theme
# 更改conf.py
import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# 2. 代码块一键复制按钮插件
pip install -U sphinx_copybutton
# 3. mermaid图表渲染插件
pip install -U sphinxcontrib.mermaid
# 4.



