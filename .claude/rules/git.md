# Git 규칙

## 작업 흐름
- `main` 브랜치에 직접 커밋한다.
- 작업 시작 전 `git pull`, 커밋 후 `git push`.
- Claude는 사용자가 요청할 때만 커밋하고, **push는 사용자 확인 후** 실행한다.
- `git push --force`, `git reset --hard` 등 기록을 지우는 명령은 쓰지 않는다. 충돌이 나면 멈추고 사용자에게 알린다.

## 커밋 메시지 (한글)
- 형식: `[분류] 무엇을 했는지` (한 줄, 명사형으로 끝냄)
- 분류: `회의록`, `주제선정`, `보고서`, `설정`, `기타`
- 예시
  - `[회의록] 2026-09-28 회의록 추가`
  - `[주제선정] 후보주제 목록에 예인선 전복 방지 장치 추가`
  - `[설정] .claude 규칙·스킬 추가`

## 올리지 않는 것
- `.env`, `.DS_Store`, `.claude/settings.local.json` (개인 설정)
- 강의노트·강의녹음·제출 PDF 원본, 개인 아이디어 자료

## Mac/Windows 혼용
- 한글 파일명이 Mac에서 자모 분리(NFD)로 올라가지 않도록 Mac에서는 `git config core.precomposeunicode true`.
- 한글 경로를 읽기 좋게 보려면 `git config core.quotepath false`.
- Windows에서는 줄바꿈 변환을 위해 `git config core.autocrlf true`, Mac에서는 `input`.
