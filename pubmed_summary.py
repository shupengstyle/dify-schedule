import os
import time
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from Bio import Entrez
from bs4 import BeautifulSoup
import google.generativeai as genai
import requests
import logging
import dotenv
import json
from datetime import datetime, timedelta

# Load environment variables from .env file
dotenv.load_dotenv()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 环境变量配置
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY")
EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")
EMAIL_PASSWORD = os.getenv("EMAIL_PASSWORD")
SMTP_SERVER = os.getenv("EMAIL_SMTP_SERVER", "smtp.yeah.net")
SMTP_PORT = int(os.getenv("EMAIL_SMTP_PORT", 465))
SEARCH_QUERY = os.getenv("SEARCH_QUERY")
MAX_RESULTS = int(os.getenv("MAX_RESULTS", 5))
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")
SUMMARY_LANGUAGE = os.getenv("SUMMARY_LANGUAGE", "en")  # 可以配置总结语言
# File to store processed PMIDs
PROCESSED_PMIDS_FILE = os.getenv("PROCESSED_PMIDS_FILE", "processed_pmids.json")
# Expiration time for processed PMIDs (in days)
PROCESSED_PMIDS_EXPIRATION_DAYS = int(os.getenv("PROCESSED_PMIDS_EXPIRATION_DAYS", 30))

# 初始化配置
Entrez.email = EMAIL_ADDRESS
Entrez.api_key = PUBMED_API_KEY

# 配置 Gemini API
genai.configure(api_key=GOOGLE_API_KEY)
model = genai.GenerativeModel('gemini-pro')


def load_processed_pmids():
    """Load processed PMIDs from file."""
    try:
        with open(PROCESSED_PMIDS_FILE, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        return []
    except json.JSONDecodeError:  # Handle corrupted JSON file
        logging.warning("Processed PMIDs file is corrupted. Starting with an empty list.")
        return []


def save_processed_pmids(pmids):
    """Save processed PMIDs to file."""
    try:
        with open(PROCESSED_PMIDS_FILE, "w") as f:
            json.dump(pmids, f)
    except Exception as e:
        logging.error(f"Failed to save processed PMIDs to file: {e}")


def fetch_articles():
    """获取文献列表"""
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=SEARCH_QUERY,
            retmax=MAX_RESULTS,
            sort="pub_date",
            usehistory="y"
        )
        result = Entrez.read(handle)
        return result["IdList"]
    except Exception as e:
        logging.error(f"文献搜索失败: {str(e)}")  # 使用logging
        return []


def translate_text(text, target_language="zh-CN"):
    """使用 Gemini API 翻译文本"""
    if not text:
        return ""
    try:
        response = model.generate_content(f"Translate the following text to {target_language}: {text}")
        return response.text.strip() if response.text else text
    except Exception as e:
        logging.error(f"Gemini API 翻译失败: {str(e)}")  # 使用logging
        return text


def summarize_text(text, target_language="en"):  # 默认英文总结
    """使用 Gemini API 生成文本总结"""
    if not text:
        return "无内容可总结"
    try:
        prompt = f"""
        请以学术角度总结以下医学研究文章，去除任何与标题、作者、摘要、方法、结论等相关的明确字眼。
        请概括文章的研究目的、采用的主要方法、取得的关键研究结果以及研究的意义和价值。
        总结应简洁明了，避免冗余信息。
        文本：
        {text}
        """
        response = model.generate_content(prompt)
        # 清理 Gemini API 返回的文本，去除多余的星号和其他无关符号
        cleaned_text = response.text.strip().replace("**", "").replace("*", "",).replace("■", "").replace("●",
                                                                                                           "").replace(
            "◆", "")
        return cleaned_text if response.text else "无法生成总结"
    except Exception as e:
        logging.error(f"Gemini API 总结失败: {str(e)}")
        return "无法生成总结"


def cleanup_processed_pmids(pmids):
    """Cleanup processed PMIDs, removing entries older than PROCESSED_PMIDS_EXPIRATION_DAYS."""
    cutoff_date = datetime.now() - timedelta(days=PROCESSED_PMIDS_EXPIRATION_DAYS)
    cleaned_pmids = []
    for entry in pmids:
        try:
            # Assuming the timestamp is stored as a string in ISO format
            entry_date = datetime.fromisoformat(entry["timestamp"])
            if entry_date >= cutoff_date:
                cleaned_pmids.append(entry)
        except (KeyError, ValueError) as e:
            logging.warning(f"Invalid entry format in processed PMIDs: {entry}. Error: {e}. Skipping.")
            continue  # Skip invalid entries
    logging.info(f"Cleaned up processed PMIDs, removed {len(pmids) - len(cleaned_pmids)} entries.")
    return cleaned_pmids


# ... 其他函数 ...

if __name__ == "__main__":
    logging.info("🚀 开始获取文献...")
    article_ids = fetch_articles()

    if not article_ids:
        logging.warning("❌ 未找到相关文献")  # 使用logging
        exit()

    processed_pmids = load_processed_pmids()
    # Clean up expired PMIDs before processing
    processed_pmids = cleanup_processed_pmids(processed_pmids)

    new_articles = []

    for pmid in article_ids:
        # Format processed_pmids correctly for checking
        processed_pmids_ids = [entry['pmid'] for entry in processed_pmids]

        if pmid not in processed_pmids_ids:
            try:
                if article := get_article_details(pmid):
                    new_articles.append(article)
                    # Store PMID with timestamp
                    processed_pmids.append({"pmid": pmid, "timestamp": datetime.now().isoformat()})
                time.sleep(0.5)  # 避免速率限制，根据实际情况调整
            except Exception as e:
                logging.exception(f"处理PMID {pmid} 时发生未知错误")  # 记录完整堆栈信息
        else:
            logging.info(f"PMID {pmid} 已经处理过，跳过")

    if new_articles:
        send_email(new_articles)  # 只发送新的文献
        save_processed_pmids(processed_pmids)  # 保存已处理的PMID列表
    else:
        logging.info("❌ 没有新的文献需要发送")
