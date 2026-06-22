# Sigma Protocol v3.0 — claude_call은 비활성화됨
# 판단/취합은 각 에이전트에서 순수 Python 로직으로 직접 처리.
# API 키 없이 Claude Code 환경에서 동작.


def call_claude(prompt: str, system: str = None, **kwargs) -> str:
    raise NotImplementedError(
        "call_claude()는 비활성화되었습니다. "
        "각 에이전트는 순수 Python 룰베이스 로직으로 판단합니다."
    )
