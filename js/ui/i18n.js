// Interface languages for navigation and common controls. Engineering content (equation names,
// analysis titles, units) stays in English, the working language of aviation documentation.

import { state } from '../core/store.js';

export const LANGS = { en: 'English', fr: 'Français', es: 'Español', pt: 'Português', de: 'Deutsch', ar: 'العربية', zh: '中文', hi: 'हिन्दी', ru: 'Русский', sw: 'Kiswahili' };
const D = {
  'Overview': { fr: 'Vue d’ensemble', es: 'Resumen', pt: 'Visão geral', de: 'Übersicht', ar: 'نظرة عامة', zh: '概览', hi: 'अवलोकन', ru: 'Обзор', sw: 'Muhtasari' },
  'Case & input portal': { fr: 'Cas et données d’entrée', es: 'Caso y datos de entrada', pt: 'Caso e dados de entrada', de: 'Fall und Eingabedaten', ar: 'الحالة وبيانات الإدخال', zh: '算例与输入', hi: 'केस और इनपुट', ru: 'Задача и исходные данные', sw: 'Kesi na data za kuingiza' },
  'Geometry & mesh': { fr: 'Géométrie et maillage', es: 'Geometría y malla', pt: 'Geometria e malha', de: 'Geometrie und Netz', ar: 'الهندسة والشبكة', zh: '几何与网格', hi: 'ज्यामिति और मेश', ru: 'Геометрия и сетка', sw: 'Jiometri na matundu' },
  'Materials library': { fr: 'Bibliothèque de matériaux', es: 'Biblioteca de materiales', pt: 'Biblioteca de materiais', de: 'Werkstoffbibliothek', ar: 'مكتبة المواد', zh: '材料库', hi: 'सामग्री पुस्तकालय', ru: 'Библиотека материалов', sw: 'Maktaba ya nyenzo' },
  'Integrated run': { fr: 'Calcul intégré', es: 'Ejecución integrada', pt: 'Execução integrada', de: 'Gekoppelter Lauf', ar: 'تشغيل متكامل', zh: '综合运行', hi: 'एकीकृत रन', ru: 'Связанный расчёт', sw: 'Uendeshaji jumuishi' },
  'High-fidelity bridge': { fr: 'Passerelle haute fidélité', es: 'Puente de alta fidelidad', pt: 'Ponte de alta fidelidade', de: 'High-Fidelity-Brücke', ar: 'جسر الدقة العالية', zh: '高保真接口', hi: 'उच्च-निष्ठा सेतु', ru: 'Мост к точным решателям', sw: 'Daraja la usahihi wa juu' },
  'Decision support': { fr: 'Aide à la décision', es: 'Apoyo a la decisión', pt: 'Apoio à decisão', de: 'Entscheidungshilfe', ar: 'دعم القرار', zh: '决策支持', hi: 'निर्णय सहायता', ru: 'Поддержка решений', sw: 'Usaidizi wa maamuzi' },
  'Live data': { fr: 'Données en direct', es: 'Datos en vivo', pt: 'Dados ao vivo', de: 'Live-Daten', ar: 'بيانات مباشرة', zh: '实时数据', hi: 'लाइव डेटा', ru: 'Живые данные', sw: 'Data za moja kwa moja' },
  'Reports & assurance': { fr: 'Rapports et assurance', es: 'Informes y garantía', pt: 'Relatórios e garantia', de: 'Berichte und Nachweise', ar: 'التقارير والضمان', zh: '报告与验证', hi: 'रिपोर्ट और आश्वासन', ru: 'Отчёты и подтверждение', sw: 'Ripoti na uhakikisho' },
  'Install, offline & about': { fr: 'Installer, hors ligne, à propos', es: 'Instalar, sin conexión y acerca de', pt: 'Instalar, offline e sobre', de: 'Installieren, offline, Info', ar: 'التثبيت وبدون اتصال وحول', zh: '安装、离线与关于', hi: 'इंस्टॉल, ऑफ़लाइन और परिचय', ru: 'Установка, офлайн, о программе', sw: 'Sakinisha, nje ya mtandao na kuhusu' },
  'Install': { fr: 'Installer', es: 'Instalar', pt: 'Instalar', de: 'Installieren', ar: 'تثبيت', zh: '安装', hi: 'इंस्टॉल', ru: 'Установить', sw: 'Sakinisha' },
  'Previous': { fr: 'Précédent', es: 'Anterior', pt: 'Anterior', de: 'Zurück', ar: 'السابق', zh: '上一页', hi: 'पिछला', ru: 'Назад', sw: 'Iliyotangulia' },
  'Next': { fr: 'Suivant', es: 'Siguiente', pt: 'Seguinte', de: 'Weiter', ar: 'التالي', zh: '下一页', hi: 'अगला', ru: 'Далее', sw: 'Inayofuata' },
  'Suites': { fr: 'Modules', es: 'Módulos', pt: 'Módulos', de: 'Module', ar: 'الوحدات', zh: '模块', hi: 'सूट', ru: 'Модули', sw: 'Moduli' },
  'Case': { fr: 'Cas', es: 'Caso', pt: 'Caso', de: 'Fall', ar: 'الحالة', zh: '算例', hi: 'केस', ru: 'Задача', sw: 'Kesi' },
  'Run all': { fr: 'Tout lancer', es: 'Ejecutar todo', pt: 'Executar tudo', de: 'Alles rechnen', ar: 'تشغيل الكل', zh: '全部运行', hi: 'सब चलाएँ', ru: 'Запустить всё', sw: 'Endesha zote' },
  'Advice': { fr: 'Conseils', es: 'Consejos', pt: 'Conselhos', de: 'Empfehlung', ar: 'توصيات', zh: '建议', hi: 'सलाह', ru: 'Советы', sw: 'Ushauri' },
  'Run analysis': { fr: 'Lancer l’analyse', es: 'Ejecutar análisis', pt: 'Executar análise', de: 'Analyse starten', ar: 'تشغيل التحليل', zh: '运行分析', hi: 'विश्लेषण चलाएँ', ru: 'Запустить расчёт', sw: 'Endesha uchambuzi' },
  'Run whole suite': { fr: 'Lancer tout le module', es: 'Ejecutar todo el módulo', pt: 'Executar todo o módulo', de: 'Ganzes Modul rechnen', ar: 'تشغيل الوحدة كاملة', zh: '运行整个模块', hi: 'पूरा सूट चलाएँ', ru: 'Запустить весь модуль', sw: 'Endesha moduli yote' },
  'Set up & run': { fr: 'Configurer et lancer', es: 'Configurar y ejecutar', pt: 'Configurar e executar', de: 'Einrichten und rechnen', ar: 'الإعداد والتشغيل', zh: '设置并运行', hi: 'सेटअप और रन', ru: 'Настройка и запуск', sw: 'Andaa na endesha' },
  'Mesh & convergence': { fr: 'Maillage et convergence', es: 'Malla y convergencia', pt: 'Malha e convergência', de: 'Netz und Konvergenz', ar: 'الشبكة والتقارب', zh: '网格与收敛', hi: 'मेश और अभिसरण', ru: 'Сетка и сходимость', sw: 'Matundu na muunganiko' },
  'Studies': { fr: 'Études', es: 'Estudios', pt: 'Estudos', de: 'Studien', ar: 'دراسات', zh: '研究', hi: 'अध्ययन', ru: 'Исследования', sw: 'Tafiti' },
  'Verification & validation': { fr: 'Vérification et validation', es: 'Verificación y validación', pt: 'Verificação e validação', de: 'Verifikation und Validierung', ar: 'التحقق والمصادقة', zh: '验证与确认', hi: 'सत्यापन और प्रमाणीकरण', ru: 'Верификация и валидация', sw: 'Uhakiki na uthibitisho' },
  'Specification': { fr: 'Spécification', es: 'Especificación', pt: 'Especificação', de: 'Spezifikation', ar: 'المواصفات', zh: '规范', hi: 'विनिर्देश', ru: 'Спецификация', sw: 'Vipimo' },
  'Live resources': { fr: 'Ressources en direct', es: 'Recursos en vivo', pt: 'Recursos ao vivo', de: 'Live-Quellen', ar: 'مصادر مباشرة', zh: '实时资源', hi: 'लाइव संसाधन', ru: 'Актуальные источники', sw: 'Rasilimali za moja kwa moja' },
  'Set up your aircraft': { fr: 'Définir votre aéronef', es: 'Configure su aeronave', pt: 'Configure a sua aeronave', de: 'Luftfahrzeug einrichten', ar: 'إعداد طائرتك', zh: '设置您的飞行器', hi: 'अपना विमान सेट करें', ru: 'Задать воздушное судно', sw: 'Andaa ndege yako' },
  'Run all 26 suites': { fr: 'Lancer les 26 modules', es: 'Ejecutar los 26 módulos', pt: 'Executar os 26 módulos', de: 'Alle 26 Module rechnen', ar: 'تشغيل الوحدات الـ26', zh: '运行全部 26 个模块', hi: 'सभी 26 सूट चलाएँ', ru: 'Запустить все 26 модулей', sw: 'Endesha moduli zote 26' },
  'One connected workspace for every aircraft engineering analysis.': { fr: 'Un espace de travail connecté pour toutes les analyses d’ingénierie aéronautique.', es: 'Un espacio de trabajo conectado para todo análisis de ingeniería aeronáutica.', pt: 'Um espaço de trabalho ligado para todas as análises de engenharia aeronáutica.', de: 'Ein vernetzter Arbeitsplatz für jede luftfahrttechnische Analyse.', ar: 'مساحة عمل مترابطة لكل تحليلات هندسة الطائرات.', zh: '面向所有飞行器工程分析的一体化工作空间。', hi: 'हर विमान इंजीनियरिंग विश्लेषण के लिए एक जुड़ा हुआ कार्यक्षेत्र।', ru: 'Единая связанная среда для всех инженерных расчётов воздушных судов.', sw: 'Eneo moja la kazi lililounganishwa kwa kila uchambuzi wa uhandisi wa ndege.' },
  'Search suites, analyses and pages…': { fr: 'Rechercher modules, analyses et pages…', es: 'Buscar módulos, análisis y páginas…', pt: 'Procurar módulos, análises e páginas…', de: 'Module, Analysen und Seiten suchen…', ar: 'ابحث في الوحدات والتحليلات والصفحات…', zh: '搜索模块、分析和页面…', hi: 'सूट, विश्लेषण और पृष्ठ खोजें…', ru: 'Поиск модулей, расчётов и страниц…', sw: 'Tafuta moduli, uchambuzi na kurasa…' },
};
export const lang = () => state.settings.lang || 'en';
/** Translate an interface string; falls back to English. */
export const t = (s) => (lang() === 'en' ? s : D[s]?.[lang()] ?? s);
export function applyLang() { document.documentElement.lang = lang(); document.documentElement.dir = lang() === 'ar' ? 'rtl' : 'ltr'; }
