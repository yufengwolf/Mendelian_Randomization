# Mendelian_Randomization

## 前提准备

1. **操作系统**: Windows
2. **环境要求**:
   - R
   - Python 3

3. **获取 Token**:
   - 访问 [OpenGWAS API](https://api.opengwas.io/profile/) 申请 Token (这个有时效性, 大概每2周需要自己手动更新, 不然就无法读取数据库)
   - 将 Token 存储在本地用户目录下的 `.Renviron` 文件中（如果没有该文件，请自行创建）

## 操作步骤

1. 打开一个 Bash 命令行界面
2. 运行以下命令启动应用程序：
   ```bash
   python app.py
3. 等待片刻，应用程序启动后，打开浏览器并输入以下地址：
```bash
  http://127.0.0.1:5000/
浏览器将自动跳转到分析界面
参考资料
具体操作可以参考这篇文章: [CSDN 博客](https://blog.csdn.net/qq994327432/article/details/142111707?spm=1001.2014.3001.5502)
